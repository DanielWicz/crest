!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2021 - 2022 Philipp Pracht
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

!> module xtb_sc
!> A module containing routines for
!> system calls to the xtb code

!=========================================================================================!
module xtb_sc
  use iso_fortran_env,only:wp => real64
  use strucrd
  use calc_type
  use iomod,only:makedir,directory_exist,remove,command
  implicit none
  !>--- private module variables and parameters
  private
  !> Files wiped from the calculation space before every xtb call.
  !> xtbinp.engrad and xtb.out MUST be in this list: both are READ back after
  !> the run, so a leftover from the previous structure in the same per-thread
  !> calcspace could be parsed as if it belonged to the current one. (Same
  !> failure mode as the %gradfile aliasing bug fixed in src/algos/parallel.f90.)
  integer,parameter :: nf = 5
  character(len=*),parameter :: xtbfiles(nf) = [&
          & 'charges      ','xtbinp.grad  ','xtbrestart   ', &
          & 'xtbinp.engrad','xtb.out      ']
  character(len=3),parameter :: xtb = 'xtb'
  character(len=10),parameter :: xyzn = 'xtbinp.xyz'
  character(len=13),parameter :: gf = 'xtbinp.engrad'

  public :: xtb_engrad

!========================================================================================!
!========================================================================================!
contains  !>--- Module routines start here
!========================================================================================!
!========================================================================================!

  subroutine xtb_engrad(mol,calc,energy,grad,iostatus)
    use iso_fortran_env,only:wp => real64
    use strucrd
    use calc_type
    use iomod,only:makedir,directory_exist,remove

    implicit none
    type(coord) :: mol
    type(calculation_settings) :: calc

    real(wp),intent(inout) :: energy
    real(wp),intent(inout) :: grad(3,mol%nat)
    integer,intent(out) :: iostatus

    integer :: i,j,k,l,ich,och,io
    logical :: ex

    iostatus = 0

    !>--- setup system call information
    !$omp critical
    call xtb_setup(mol,calc)
    !$omp end critical

    !>--- do the systemcall
    call initsignal()
    call command(calc%systemcall,iostatus)
    if (iostatus /= 0) return

    !>--- read energy (and gradient, unless this is an energy-only job)
    !$omp critical
    if (.not.calc%energyonly) then
      call rd_xtb_engrad(mol,calc,energy,grad,iostatus)
    else
      call rd_xtb_energy(mol,calc,energy,grad,iostatus)
    end if
    !$omp end critical
    if (iostatus /= 0) return

    !>--- read WBOs?
    !$omp critical
    call rd_xtb_wbo(mol,calc,iostatus)
    !$omp end critical
    if (iostatus /= 0) return

    return
  end subroutine xtb_engrad

!========================================================================================!
  subroutine xtb_setup(mol,calc)
    use iso_fortran_env,only:wp => real64
    use strucrd
    use calc_type
    use iomod,only:makedir,directory_exist,remove

    implicit none
    type(coord) :: mol
    type(calculation_settings) :: calc

    integer :: l
    character(len=:),allocatable :: fname
    character(len=:),allocatable :: cpath
    character(len=10) :: num
    integer :: i,j,k,ich,och,io
    logical :: ex

    call initsignal()

    !>--- set default binary if not present
    if (.not.allocated(calc%binary)) then
      calc%binary = xtb
    end if

    !>--- check for the calculation space
    if (allocated(calc%calcspace)) then
      ex = directory_exist(calc%calcspace)
      if (.not.ex) then
        io = makedir(trim(calc%calcspace))
      end if
      cpath = calc%calcspace
    else
      cpath = ''
    end if
    !>--- cleanup old files
    do i = 1,nf
      !write(*,*) trim(cpath)//sep//trim(xtbfiles(i))
      call remove(trim(cpath)//sep//trim(xtbfiles(i)))
    end do
    deallocate (cpath)

    !>--- construct path information and write coord file
    if (.not.allocated(calc%calcfile)) then
      if (allocated(calc%calcspace)) then
        l = len_trim(calc%calcspace)
        fname = trim(calc%calcspace)
        if (calc%calcspace(l:l) == sep) then
          fname = trim(fname)//xyzn
        else
          fname = trim(fname)//sep//xyzn
        end if
      else
        fname = xyzn
      end if
      calc%calcfile = fname
    else
      fname = calc%calcfile
    end if
    call mol%write(fname)
    deallocate (fname)

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
    !>--- if the systemcall was already set up, return
    if (allocated(calc%systemcall)) return
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!

    !>--- construct path information for sys-call
    if (allocated(calc%calcspace)) then
      calc%systemcall = 'cd '//calc%calcspace//' &&'
      calc%systemcall = trim(calc%systemcall)//' '//trim(calc%binary)
    else
      calc%systemcall = trim(calc%binary)
    end if

    !>--- add other call information
    calc%systemcall = trim(calc%systemcall)//' '//xyzn
    !>--- chrg and uhf
    if (calc%chrg .ne. 0) then
      write (num,'(i0)') calc%chrg
      calc%systemcall = trim(calc%systemcall)//' '//'--chrg'
      calc%systemcall = trim(calc%systemcall)//' '//trim(num)
    end if
    if (calc%uhf .ne. 0) then
      write (num,'(i0)') calc%uhf
      calc%systemcall = trim(calc%systemcall)//' '//'--uhf'
      calc%systemcall = trim(calc%systemcall)//' '//trim(num)
    end if
    !>--- user-set flags
    if (allocated(calc%other)) then
      calc%systemcall = trim(calc%systemcall)//' '//trim(calc%other)
    end if
    !>--- implicit solvation
    if (allocated(calc%solvmodel).and.allocated(calc%solvent)) then
      select case (trim(calc%solvmodel))
      case ('gbsa','alpb','cpcm')
        if (index(calc%systemcall,'--'//trim(calc%solvmodel)) .eq. 0) then
          calc%systemcall = trim(calc%systemcall)//' --'//trim(calc%solvmodel)
          calc%systemcall = trim(calc%systemcall)//' '//trim(calc%solvent)
        end if
      end select
      !>--- DRACO charge-dependent cavity radii. Only meaningful together with an
      !>--- implicit solvation model, so it is appended inside this branch.
      if (allocated(calc%draco)) then
        if (index(calc%systemcall,'--draco') .eq. 0) then
          calc%systemcall = trim(calc%systemcall)//' --draco'
          if (len_trim(calc%draco) > 0) then
            calc%systemcall = trim(calc%systemcall)//' '//trim(calc%draco)
          end if
        end if
      end if
    end if
    !>--- gradient or energy only?
    !> The refinement singlepoint loop (crest_sploop) throws the gradient away,
    !> so calculators flagged rdgrad=.false. ask xtb not to compute one. This is
    !> free for GFN1/GFN2, whose analytic gradient costs almost nothing on top
    !> of the SCF, but g-xTB computes an analytic CPHF/Z-vector gradient on
    !> every singlepoint: measured on 50 atoms, 0.47 s with --grad vs 0.30 s
    !> with --sp-nograd, same energy to all 12 printed digits.
    !> NOTE on the duplicate guards: '--sp-nograd' does NOT contain the
    !> substring '-grad' (the character before "grad" is 'o'), so the two
    !> branches need different tests. Checking for '-grad' in the energy-only
    !> branch would be dead code.
    if (.not.calc%energyonly) then
      !>--- don't miss the --grad flag!
      if (index(calc%systemcall,'-grad') .eq. 0) then
        calc%systemcall = trim(calc%systemcall)//' '//'--grad'
      end if
    else
      if (index(calc%systemcall,'nograd') .eq. 0) then
        calc%systemcall = trim(calc%systemcall)//' '//'--sp-nograd'
      end if
    end if

    !>--- add printout information
    calc%systemcall = trim(calc%systemcall)//' '//'> xtb.out'
    calc%systemcall = trim(calc%systemcall)//dev0

    !write (*,*) calc%systemcall
    return
  end subroutine xtb_setup

!========================================================================================!
! subroutine rd_xtb_energy
! Read ONLY the total energy, from xtb's stdout log.
!
! Used for energy-only jobs (calc%energyonly), where xtb was called with
! --sp-nograd and therefore wrote no .engrad file. The printed
!     | TOTAL ENERGY   <E> Eh |
! line carries the same 12 decimals as the .engrad file (verified identical on
! GFN2 and g-xTB), which is far more than the ~1e-6 Eh that ensemble ranking
! needs.
!
! Safety rules:
!   1. the FIRST TOTAL ENERGY match wins and the scan stops there. xtb prints
!      exactly one such line per run (property.F90 write_energy /
!      write_energy_gff), so this is complete, and it keeps the read - which
!      sits inside a global !$omp critical - from walking a large g-xTB log;
!   2. xtb.out is wiped by xtb_setup before every call (see xtbfiles), so a
!      stale log from the previous structure in a reused per-thread calcspace
!      cannot be mistaken for this one;
!   3. the caller only reaches this routine after the xtb systemcall returned 0
!      (see xtb_engrad), which is the same guarantee rd_xtb_engrad relies on.
! NOTE: do NOT try to also require xtb's "normal termination of xtb" banner.
! xtb prints it on STDERR, and the systemcall sends stderr to /dev/null (the
! dev0 parameter), so it is never present in xtb.out and the check would fail
! every single time.
  subroutine rd_xtb_energy(mol,calc,energy,grad,iostatus)
    use iso_fortran_env,only:wp => real64
    use ieee_arithmetic,only:ieee_is_finite
    use strucrd
    use calc_type
    implicit none
    type(coord) :: mol
    type(calculation_settings) :: calc
    real(wp),intent(inout) :: energy
    real(wp),intent(inout) :: grad(3,mol%nat)
    integer,intent(out) :: iostatus

    character(len=:),allocatable :: outfile
    character(len=256) :: atmp
    integer :: ich,io,k
    logical :: ex,found
    real(wp) :: edum

    call initsignal()
    iostatus = 0
    energy = 0.0_wp
    grad = 0.0_wp

    if (allocated(calc%calcspace)) then
      outfile = trim(calc%calcspace)//sep//'xtb.out'
    else
      outfile = 'xtb.out'
    end if

    inquire (file=outfile,exist=ex)
    if (.not.ex) then
      iostatus = 1
      return
    end if

    found = .false.
    open (newunit=ich,file=outfile,status='old',action='read',iostat=io)
    if (io /= 0) then
      iostatus = 1
      return
    end if
    do
      read (ich,'(a)',iostat=io) atmp
      if (io /= 0) exit
      k = index(atmp,'TOTAL ENERGY')
      if (k > 0) then
        !> layout: "| TOTAL ENERGY  <value> Eh   |"
        read (atmp(k+12:),*,iostat=io) edum
        if (io == 0) then
          !> A NaN energy passes list-directed input happily, and xtb can print
          !> "TOTAL ENERGY NaN Eh" while exiting 0. Reject it here; the .engrad
          !> route has no equivalent trap either, but this one is new code.
          if (ieee_is_finite(edum)) then
            energy = edum
            found = .true.
            exit
          end if
        end if
      end if
    end do
    close (ich)

    if (.not.found) then
      energy = 0.0_wp
      iostatus = 1
    end if

    return
  end subroutine rd_xtb_energy

!========================================================================================!
! subroutine rd_xtb_engrad
! read xtb's energy and Cartesian gradient from file
! xtb's *.engrad format is used for this
  subroutine rd_xtb_engrad(mol,calc,energy,grad,iostatus)
    use iso_fortran_env,only:wp => real64
    use strucrd
    use calc_type
    use iomod,only:makedir,directory_exist,remove
    use gradreader_module,only:rd_grad_engrad

    implicit none
    type(coord) :: mol
    type(calculation_settings) :: calc
    real(wp),intent(inout) :: energy
    real(wp),intent(inout) :: grad(3,mol%nat)
    integer,intent(out) :: iostatus
    integer :: n,c
    real(wp) :: dum
    character(len=128) :: atmp

    integer :: i,j,k,ich,och,io
    logical :: ex

    call initsignal()

    iostatus = 0

    if (.not.allocated(calc%gradfile)) then
      if (allocated(calc%calcspace)) then
        calc%gradfile = trim(calc%calcspace)//sep//gf
      else
        calc%gradfile = gf
      end if
    end if

    inquire (file=calc%gradfile,exist=ex)
    if (.not.ex) then
      iostatus = 1
      return
    end if

    c = 0
    open (newunit=ich,file=calc%gradfile)
    call rd_grad_engrad(ich,mol%nat,energy,grad,iostatus)
    close (ich)

    return
  end subroutine rd_xtb_engrad
!========================================================================================!
! subroutine rd_xtb_wbo
! helper routine ro read xtb WBOs
  subroutine rd_xtb_wbo(mol,calc,iostatus)
    implicit none
    type(coord) :: mol
    type(calculation_settings) :: calc
    integer,intent(out) :: iostatus

    real(wp) :: dum
    character(len=:),allocatable :: wbofile
    character(len=128) :: atmp

    integer :: i,j,k,l,ich,och,io
    logical :: ex
    call initsignal()

    iostatus = 0

    if (calc%rdwbo) then
      if (allocated(calc%calcspace)) then
        wbofile = trim(calc%calcspace)//sep//'wbo'
      else
        wbofile = 'wbo'
      end if
    else
      return
    end if

    inquire (file=wbofile,exist=ex)
    if (.not.ex) then
      iostatus = 1
      return
    end if

    if (allocated(calc%wbo)) deallocate (calc%wbo)
    allocate (calc%wbo(mol%nat,mol%nat),source=0.0_wp)

    open (newunit=ich,file=wbofile)
    do
      read (ich,'(a)',iostat=io) atmp
      if (io < 0) exit
      read (atmp,*) i,j,dum
      calc%wbo(i,j) = dum
      calc%wbo(j,i) = dum
    end do
    close (ich)

  end subroutine rd_xtb_wbo

!========================================================================================!
end module xtb_sc
