!================================================================================!
! This file is part of crest.
!
! crest is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! crest is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with crest.  If not, see <https://www.gnu.org/licenses/>.
!================================================================================!

!> Gas-phase g-xTB re-optimiser used by "-gfnff:gxtb".
!>
!> One native "xtb --gxtb --opt <level>" process is run per structure, at most
!> `njobs` of them at a time, every one single-threaded and unconstrained
!> (no CREST calculator, no restraint, no wall potential, no solvent).
!> Results are stored by structure index, never in completion order, so the
!> outcome does not depend on the number of jobs or on scheduling.
!>
!> Encapsulation: all state of the optimiser is private. The caller sees
!> init / reoptimize and read-only queries (nfailed, failed_indices).
module gxtb_reopt_module
  use crest_parameters,only:wp,stdout
  use omp_lib
  use strucrd,only:i2e
  implicit none
  private

  public :: gxtb_reoptimizer
  public :: gxtb_level_is_valid

  !> Energy assigned to a structure whose optimisation failed (CREST's
  !> failure convention: any positive value lies above every real total energy).
  real(wp),parameter,public :: gxtb_efail = 1.0_wp

  integer,parameter :: nlevels = 8
  character(len=6),parameter :: level_words(nlevels) = &
  & [character(len=6) :: 'crude','sloppy','loose','lax','normal','tight','vtight','extreme']

  type :: gxtb_reoptimizer
    private
    character(len=:),allocatable :: binary    !> absolute path of the xtb executable
    character(len=:),allocatable :: level     !> xtb optimisation level word
    character(len=:),allocatable :: scratch   !> scratch root directory
    integer :: charge = 0
    integer :: uhf = 0
    integer :: njobs = 1                      !> parallel single-thread jobs
    logical :: keep = .false.                 !> keep scratch directories
    logical :: ready = .false.
    integer,allocatable :: failed(:)          !> 1-based indices of failed structures
  contains
    procedure,public :: init => reopt_init
    procedure,public :: reoptimize => reopt_run
    procedure,public :: nfailed => reopt_nfailed
    procedure,public :: failed_indices => reopt_failed_indices
    procedure,private :: run_one => reopt_run_one
  end type gxtb_reoptimizer

contains
!================================================================================!

  pure function gxtb_level_is_valid(word) result(valid)
    character(len=*),intent(in) :: word
    logical :: valid
    valid = any(level_words == trim(word))
  end function gxtb_level_is_valid

!================================================================================!

  subroutine reopt_init(self,binary,charge,uhf,level,njobs,scratch,keep,ok)
    !> Configure the re-optimiser. ok = .false. if a setting is unusable.
    class(gxtb_reoptimizer),intent(inout) :: self
    character(len=*),intent(in) :: binary
    integer,intent(in) :: charge,uhf
    character(len=*),intent(in) :: level
    integer,intent(in) :: njobs
    character(len=*),intent(in) :: scratch
    logical,intent(in) :: keep
    logical,intent(out) :: ok
    ok = .false.
    self%ready = .false.
    if (len_trim(binary) < 1 .or. len_trim(scratch) < 1) return
    if (.not.gxtb_level_is_valid(level)) return
    if (njobs < 1 .or. uhf < 0) return
    self%binary = trim(binary)
    self%level = trim(level)
    self%scratch = trim(scratch)
    self%charge = charge
    self%uhf = uhf
    self%njobs = njobs
    self%keep = keep
    if (allocated(self%failed)) deallocate (self%failed)
    allocate (self%failed(0))
    self%ready = .true.
    ok = .true.
  end subroutine reopt_init

!================================================================================!

  subroutine reopt_run(self,nat,nall,at,xyz_ang,energies,ok)
    !> Re-optimise all structures. xyz_ang (Angstrom) is overwritten by the
    !> optimised geometry of every structure that converged and is left as
    !> given for the others, which get energies = gxtb_efail and ok = .false.
    class(gxtb_reoptimizer),intent(inout) :: self
    integer,intent(in) :: nat,nall
    integer,intent(in) :: at(nat)
    real(wp),intent(inout) :: xyz_ang(3,nat,nall)
    real(wp),intent(out) :: energies(nall)
    logical,intent(out) :: ok(nall)
    integer :: k,nj,io
    logical :: good
    integer :: nfail
    ok = .false.
    energies = gxtb_efail
    if (allocated(self%failed)) deallocate (self%failed)
    allocate (self%failed(0))
    if (.not.self%ready .or. nall < 1) return
    call execute_command_line("mkdir -p '"//self%scratch//"'",exitstat=io)
    if (io /= 0) return
    nj = max(1,min(self%njobs,nall))
!$omp parallel do num_threads(nj) schedule(dynamic,1) default(shared) private(k,good)
    do k = 1,nall
      call self%run_one(k,nat,at,xyz_ang(:,:,k),energies(k),good)
      ok(k) = good
    end do
!$omp end parallel do
    nfail = count(.not.ok)
    deallocate (self%failed)
    allocate (self%failed(nfail))
    nfail = 0
    do k = 1,nall   !> ascending order, independent of completion order
      if (.not.ok(k)) then
        nfail = nfail+1
        self%failed(nfail) = k
      end if
    end do
    if (.not.self%keep) then
      call execute_command_line("rm -rf '"//self%scratch//"'",exitstat=io)
    end if
  end subroutine reopt_run

!================================================================================!

  subroutine reopt_run_one(self,k,nat,at,xyz,energy,good)
    !> One structure, own directory <scratch>/s<k>. xyz is updated on success.
    class(gxtb_reoptimizer),intent(in) :: self
    integer,intent(in) :: k,nat
    integer,intent(in) :: at(nat)
    real(wp),intent(inout) :: xyz(3,nat)
    real(wp),intent(out) :: energy
    logical,intent(out) :: good
    character(len=:),allocatable :: dir,cmd
    character(len=16) :: ctmp
    character(len=32) :: qstr,ustr
    integer :: io,stat,ich,j
    logical :: ex
    real(wp) :: xnew(3,nat),enew
    energy = gxtb_efail
    good = .false.
    write (ctmp,'(i0)') k
    dir = self%scratch//'/s'//trim(ctmp)
    call execute_command_line("rm -rf '"//dir//"' && mkdir -p '"//dir//"'",exitstat=io)
    if (io /= 0) return
    !> input structure, Angstrom, full CREST precision
    open (newunit=ich,file=dir//'/xtbinp.xyz',status='replace',iostat=io)
    if (io /= 0) return
    write (ich,'(2x,i0)') nat
    write (ich,'(a)') ' input'
    do j = 1,nat
      write (ich,'(1x,a2,1x,3f20.10)') i2e(at(j),'nc'),xyz(1:3,j)
    end do
    close (ich)
    write (qstr,'(i0)') self%charge
    write (ustr,'(i0)') self%uhf
    cmd = "cd '"//dir//"' && OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 '"// &
    & self%binary//"' xtbinp.xyz --gxtb --chrg "//trim(qstr)//" --uhf "//trim(ustr)// &
    & " --opt "//self%level//" > xtb.out 2> xtb.err"
    call execute_command_line(cmd,exitstat=stat)
    if (stat /= 0) return
    inquire (file=dir//'/.xtboptok',exist=ex)
    if (.not.ex) return
    !> xtb 6.7.1 writes .xtboptok and exits 0 even when it stops at the cycle
    !> limit ("FAILED TO CONVERGE GEOMETRY OPTIMIZATION"), so the log decides.
    if (.not.file_contains(dir//'/xtb.out','GEOMETRY OPTIMIZATION CONVERGED')) return
    call read_xtbopt(dir//'/xtbopt.xyz',nat,xnew,enew,ex)
    if (.not.ex) return
    xyz = xnew
    energy = enew
    good = .true.
  end subroutine reopt_run_one

!================================================================================!

  function file_contains(fname,text) result(found)
    character(len=*),intent(in) :: fname,text
    logical :: found
    integer :: ich,io
    character(len=512) :: line
    found = .false.
    open (newunit=ich,file=fname,status='old',iostat=io)
    if (io /= 0) return
    do
      read (ich,'(a)',iostat=io) line
      if (io /= 0) exit
      if (index(line,text) > 0) then
        found = .true.
        exit
      end if
    end do
    close (ich)
  end function file_contains

  subroutine read_xtbopt(fname,nat,xyz,energy,ok)
    !> parse xtbopt.xyz: geometry (Angstrom) and "energy: <E> gnorm: <g>" comment
    character(len=*),intent(in) :: fname
    integer,intent(in) :: nat
    real(wp),intent(out) :: xyz(3,nat)
    real(wp),intent(out) :: energy
    logical,intent(out) :: ok
    integer :: ich,io,n,j,p
    character(len=512) :: line
    character(len=8) :: sym
    ok = .false.
    energy = gxtb_efail
    xyz = 0.0_wp
    open (newunit=ich,file=fname,status='old',iostat=io)
    if (io /= 0) return
    read (ich,*,iostat=io) n
    if (io /= 0 .or. n /= nat) then
      close (ich)
      return
    end if
    read (ich,'(a)',iostat=io) line
    if (io /= 0) then
      close (ich)
      return
    end if
    p = index(line,'energy:')
    if (p == 0) then
      close (ich)
      return
    end if
    read (line(p+7:),*,iostat=io) energy
    if (io /= 0) then
      close (ich)
      return
    end if
    do j = 1,nat
      read (ich,*,iostat=io) sym,xyz(1,j),xyz(2,j),xyz(3,j)
      if (io /= 0) then
        close (ich)
        return
      end if
    end do
    close (ich)
    !> reject NaN / non-physical output
    if (energy /= energy .or. any(xyz /= xyz) .or. energy >= 0.0_wp) return
    ok = .true.
  end subroutine read_xtbopt

!================================================================================!

  pure function reopt_nfailed(self) result(n)
    class(gxtb_reoptimizer),intent(in) :: self
    integer :: n
    n = 0
    if (allocated(self%failed)) n = size(self%failed)
  end function reopt_nfailed

  function reopt_failed_indices(self) result(idx)
    !> copy of the 1-based indices of failed structures (ascending)
    class(gxtb_reoptimizer),intent(in) :: self
    integer,allocatable :: idx(:)
    if (allocated(self%failed)) then
      idx = self%failed
    else
      allocate (idx(0))
    end if
  end function reopt_failed_indices

end module gxtb_reopt_module
