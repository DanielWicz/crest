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

!> Post-sampling stage of "-gfnff:gxtb" and the "--gxtbreopt" replay runtype.
!>
!> The stage takes the final GFN-FF-optimised ensemble, keeps its unique
!> conformers (CREGEN at GFN-FF energies), re-optimises them in the gas phase
!> with g-xTB (see gxtb_reopt_module) and writes the result back with g-xTB
!> energies, ready for the final CREGEN. The full search and the replay use
!> the very same routine, so they cannot drift apart.
module gfnff_gxtb_module
  use crest_parameters,only:wp,stdout,stderr
  use crest_data
  use strucrd
  use iomod
  use utilities
  use cregen_interface
  use gxtb_reopt_module
  implicit none
  private

  public :: crest_gfnff_gxtb_stage
  public :: crest_gxtb_replay
  public :: gfnff_gxtb_stage_count

  !> the stored input of the post-sampling stage (GFN-FF energies)
  character(len=*),parameter :: stored_ensemble = 'crest_gfnff_ensemble.xyz'

  !> how often the stage was entered in this program run (queried by crest_main,
  !> which insists on exactly one run for every -gfnff:gxtb workflow)
  integer,save :: nstage_runs = 0

contains
!================================================================================!

  function gfnff_gxtb_stage_count() result(n)
    integer :: n
    n = nstage_runs
  end function gfnff_gxtb_stage_count

!================================================================================!

  subroutine crest_gfnff_gxtb_stage(env,ensemble_file)
    !> Re-optimise the ensemble in "ensemble_file" with g-xTB (gas phase) and
    !> overwrite it with the result. A final CREGEN is the caller's job.
    type(systemdata),intent(inout) :: env
    character(len=*),intent(in) :: ensemble_file
    type(gxtb_reoptimizer) :: reopt
    integer :: nat,nall,nin,io,k
    integer,allocatable :: at(:)
    real(wp),allocatable :: xyz(:,:,:),eread(:)
    logical,allocatable :: ok(:)
    logical :: good
    integer :: t0,t1,rate
    character(len=:),allocatable :: cwd,line
    character(len=512) :: cdir
    character(len=64) :: comment
    integer,allocatable :: idx(:)

    nstage_runs = nstage_runs+1
    !>--- 1. keep the GFN-FF ensemble as the stored input of this stage
    call execute_command_line("cp -f '"//ensemble_file//"' '"//stored_ensemble//"'",exitstat=io)
    if (io /= 0) then
      write (stderr,'(a)') '**ERROR** could not store '//stored_ensemble
      env%iostatus_meta = status_failed
      return
    end if
    call rdensembleparam(ensemble_file,nat,nin)

    !>--- 2. presort: CREGEN at GFN-FF energies, unique conformers only
    call smallhead('g-xTB gas-phase re-optimization (GFN-FF sampling)')
    call newcregen(env,0,ensemble_file)
    call rename('crest_ensemble.xyz',ensemble_file)
    call remove(ensemble_file//'.sorted')
    call rdensembleparam(ensemble_file,nat,nall)
    if (nall < 1) then
      write (stderr,'(a)') '**ERROR** empty ensemble after the GFN-FF presort.'
      env%iostatus_meta = status_failed
      return
    end if
    allocate (at(nat),xyz(3,nat,nall),eread(nall),ok(nall))
    call rdensemble(ensemble_file,nat,nall,at,xyz,eread)

    !>--- 3. g-xTB re-optimization, T parallel single-thread jobs
    call getcwd(cdir)
    cwd = trim(cdir)
    call reopt%init(trim(env%ProgName),env%chrg,env%uhf,env%gxtb_optlevel,max(1,env%threads), &
    &    cwd//'/.gxtb_reopt',env%keepModef,good)
    if (.not.good) then
      write (stderr,'(a)') '**ERROR** invalid g-xTB re-optimization settings.'
      env%iostatus_meta = status_failed
      return
    end if
    call system_clock(t0,rate)
    call reopt%reoptimize(nat,nall,at,xyz,eread,ok)
    call system_clock(t1)

    !>--- 4. summary
    idx = reopt%failed_indices()
    write (stdout,'(/,1x,a)') 'g-xTB gas-phase re-optimization (xtb --gxtb --opt '//env%gxtb_optlevel//')'
    write (stdout,'(3x,a,i0)') 'structures in the GFN-FF ensemble : ',nin
    write (stdout,'(3x,a,i0)') 'unique conformers re-optimized    : ',nall
    write (stdout,'(3x,a,i0)') 'failed (energy +1.0, dropped)     : ',reopt%nfailed()
    if (size(idx) > 0) then
      line = ''
      do k = 1,size(idx)
        write (comment,'(i0)') idx(k)
        line = line//' '//trim(comment)
      end do
      write (stdout,'(3x,a)') 'failed structure indices (1-based)  :'//line
    end if
    write (stdout,'(3x,a,f10.2,a)') 'wall time                         : ', &
    &    real(t1-t0,wp)/real(rate,wp),' s'
    write (stdout,'(3x,a,i0)') 'parallel single-thread xtb jobs   : ',min(max(1,env%threads),nall)
    write (stdout,*)
    if (all(.not.ok)) then
      write (stderr,'(a)') '**ERROR** every g-xTB re-optimization failed.'
      env%iostatus_meta = status_failed
      return
    end if

    !>--- 5. write the ensemble back (g-xTB gas-phase energies)
    open (newunit=io,file=ensemble_file,status='replace')
    do k = 1,nall
      if (ok(k)) then
        write (comment,'(f22.12)') eread(k)
      else
        write (comment,'(f22.12,a)') gxtb_efail,'  !failed'
      end if
      call wrxyz(io,nat,at,xyz(:,:,k),trim(comment))
    end do
    close (io)
  end subroutine crest_gfnff_gxtb_stage

!================================================================================!

  subroutine crest_gxtb_replay(env,tim)
    !> Runtype "--gxtbreopt <crest_gfnff_ensemble.xyz>": steps 2-5 of the stage
    !> and the final CREGEN, on a working copy of the given file. The topology
    !> reference is the input structure (file "coord").
    type(systemdata),intent(inout) :: env
    type(timer),intent(inout) :: tim
    character(len=128) :: inpnam,outnam
    integer :: nat,nall,io
    logical :: ex

    call tim%start(14,'g-xTB re-optimization')
    inquire (file=env%ensemblename,exist=ex)
    if (.not.ex) then
      write (stderr,'(a)') '**ERROR** no ensemble file for --gxtbreopt: '//trim(env%ensemblename)
      env%iostatus_meta = status_input
      return
    end if
    call rdensembleparam(env%ensemblename,nat,nall)
    if (nall < 1 .or. nat /= env%nat) then
      write (stderr,'(a)') '**ERROR** the ensemble does not match the input structure (atoms) or is empty.'
      env%iostatus_meta = status_input
      return
    end if
    write (stdout,*)
    write (stdout,'(3x,a)') '====================================================='
    write (stdout,'(3x,a)') '|  g-xTB gas-phase re-optimization (replay)          |'
    write (stdout,'(3x,a)') '====================================================='
    write (stdout,'(1x,a,a)') 'Ensemble: ',trim(env%ensemblename)

    !>--- clean start: previous results must not influence this run
    call rmcres()
    call checkname_xyz(crefile,inpnam,outnam)
    call execute_command_line("cp -f '"//trim(env%ensemblename)//"' '"//trim(outnam)//"'",exitstat=io)
    if (io /= 0) then
      env%iostatus_meta = status_failed
      return
    end if
    call crest_gfnff_gxtb_stage(env,trim(outnam))
    if (env%iostatus_meta /= 0) return
    call sort_and_check(env,trim(outnam))
    write (stdout,'(/)')
    call smallhead('Final Ensemble Information')
    call V2terminating()
    call tim%stop(14)
  end subroutine crest_gxtb_replay

end module gfnff_gxtb_module
