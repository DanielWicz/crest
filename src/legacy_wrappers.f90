!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2023 Philipp Pracht
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

!================================================================================!
subroutine env2calc(env,calc,molin)
!******************************************
!* This piece of code generates a calcdata
!* object from the global settings in env
!******************************************
  use crest_parameters
  use crest_data
  use crest_calculator
  use strucrd
  use iomod
  implicit none
  !> INPUT
  type(systemdata),intent(inout) :: env
  type(coord),intent(in),optional :: molin
  !> OUTPUT
  type(calcdata) :: calc
  !> LOCAL
  type(calculation_settings) :: cal,cal2
  type(coord) :: mol

!>--- Calculator level

!>-- defaults to whatever env has selected or gfn0
  call cal%create(trim(env%gfnver))
  if (present(molin)) then
    mol = molin
  end if

  cal%uhf = env%uhf
  cal%chrg = env%chrg
  if (cal%id == jobtype%xtbsys) then
    cal%binary = trim(env%ProgName)
  end if
!>-- obtain WBOs OFF by default
  cal%rdwbo = .false.
  cal%rddip = .false.
  !> except for SP runtype (from command line!)
  if (env%crestver == crest_sp) then
    cal%rdgrad = env%gradsp
    if (cal%id .ne. jobtype%turbomole) then
      cal%rdwbo = .true.
      cal%rddip = .true.
      cal%rdqat = .true.
    else
      if (.not.env%gradsp) then
        cal%other = ''
      end if
    end if
  end if

  !> implicit solvation
  if (env%gbsa) then
    if (index(env%solv,'gbsa') .ne. 0) then
      cal%solvmodel = 'gbsa'
    else if (index(env%solv,'alpb') .ne. 0) then
      cal%solvmodel = 'alpb'
    else
      cal%solvmodel = 'unknown'
    end if
    cal%solvent = trim(env%solvent)
    !> DRACO is applied by the xtb binary on top of the solvation model
    if (allocated(env%draco)) cal%draco = env%draco
  end if

  !> do not reset parameters between calculations (opt for speed)
  cal%apiclean = .false.

  call cal%autocomplete(1)
  call calc%add(cal)

!>--- Refinement level
  if (trim(env%gfnver2) .ne. '') then
    env%gfnver2 = lowercase(env%gfnver2)
    call cal2%create(trim(env%gfnver2))

    cal2%chrg = cal%chrg
    cal2%uhf = cal%uhf
    if (cal2%id == jobtype%xtbsys) then
      cal2%binary = trim(env%ProgName)
    end if
    if (env%gbsa) then
      cal2%solvmodel = cal%solvmodel
      cal2%solvent = cal%solvent
      if (allocated(cal%draco)) cal2%draco = cal%draco
    end if

    call cal2%autocomplete(2)

    !> Which refinement the user asked for: -rsp/-refine -> singlepoint,
    !> -ropt -> a real re-optimization at gfnver2. This used to be hardwired to
    !> refine%singlepoint, which made -ropt a silent alias of -rsp.
    cal2%refine_lvl = env%refine_lvl_cli
    !> The rescoring loop (crest_sploop) discards the gradient, so do not pay
    !> for it -- but only where xtb actually honours the request. In xtb,
    !> --sp-nograd is acted on ONLY for g-xTB and only for a plain SCC run
    !> (src/prog/main.F90: `tblite%method == "gxtb" .and. set%sp_nograd .and.
    !> set%runtyp == p_run_scc`, and src/tblite/calculator.F90). For GFN1/GFN2
    !> it merely suppresses --grad, whose analytic gradient is almost free, so
    !> enabling it there would buy nothing and add a stdout-parsing path for no
    !> reason. Measured, 50 atoms x 24 structures, 8 threads, g-xTB rescoring:
    !> 1.310 s -> 0.547 s (2.4x), energies identical to all printed digits.
    !> A geoopt refinement obviously still needs gradients.
    if (allocated(cal2%other)) then
      if (index(cal2%other,'--gxtb') /= 0 .and. &
      &  (cal2%refine_lvl == refine%singlepoint .or. &
      &   cal2%refine_lvl == refine%correction)) then
        cal2%energyonly = .true.
      end if
    end if
    call calc%add(cal2)
    if (allocated(env%refine_queue)) deallocate (env%refine_queue)
    call env%addrefine(env%refine_lvl_cli)
  end if

  return
end subroutine env2calc

subroutine env2calc_setup(env)
!***********************************
!* Setup the calc object within env
!* (wrapper to get the mol object)
!***********************************
  use crest_data
  use crest_calculator
  use strucrd
  implicit none
  !> INOUT
  type(systemdata),intent(inout) :: env
  !> LOCAL
  type(calcdata) :: calc
  type(coord) :: mol
  interface
    subroutine env2calc(env,calc,molin)
      use crest_parameters
      use crest_data
      use crest_calculator
      use strucrd
      implicit none
      type(systemdata),intent(inout) :: env
      type(coord),intent(in),optional :: molin
      type(calcdata) :: calc
    end subroutine env2calc

  end interface

  call env%ref%to(mol)

  call env2calc(env,env%calc,mol)

  ! env%calc = calc
end subroutine env2calc_setup

!================================================================================!
subroutine env_apply_calclevel(env)
!*********************************************************************************
!* Re-apply env%gfnver as the METHOD (calculation level) of the already existing
!* env%calc, in place.
!*
!* WHY THIS EXISTS (fail-open guard):
!*   env%calc is built exactly once, in confparse (env2calc_setup), from the
!*   env%gfnver that was active at the end of argument parsing. Several runtypes
!*   -- QCG above all -- switch methods at RUNTIME by re-assigning env%gfnver
!*   (e.g. env%gfnver = env%ensemble_opt for -enslvl, or env%freqver for
!*   -freqlvl, see src/qcg/solvtool.f90). That assignment is only seen by the
!*   legacy (CREST <3.0) routines. On the new calculator route the method
!*   silently stayed whatever confparse had built, so CREST printed
!*   "Method for ensemble search: <X>" while actually running <Y>.
!*   This routine closes that gap. Keep it called next to EVERY runtime
!*   env%gfnver re-assignment.
!*
!* Only the level-defining fields (id/other/binary/description) are replaced.
!* Charge, spin, solvation, scratch directory, weights and all calcdata-level
!* state (constraints, wall potentials, optimizer settings) are preserved.
!*********************************************************************************
  use crest_parameters
  use crest_data
  use crest_calculator
  implicit none
  !> INOUT
  type(systemdata),intent(inout) :: env
  !> LOCAL
  integer :: j
  integer :: s_chrg,s_uhf,s_refine_lvl,s_id
  integer :: io_lvl
  real(wp) :: s_weight
  logical :: s_active,s_rdwbo,s_rddip,s_rdqat,s_rdgrad,s_apiclean
  character(len=:),allocatable :: s_calcspace,s_solvmodel,s_solvent,s_draco
  character(len=:),allocatable :: s_other,s_binary

  if (env%legacy) return
  if (env%calc%ncalculations < 1) return
  if (len_trim(env%gfnver) < 1) return

  do j = 1,env%calc%ncalculations
    !>--- save everything that create() would wipe but that is NOT part of
    !>--- the method definition itself
    s_chrg = env%calc%calcs(j)%chrg
    s_uhf = env%calc%calcs(j)%uhf
    s_refine_lvl = env%calc%calcs(j)%refine_lvl
    s_weight = env%calc%calcs(j)%weight
    s_active = env%calc%calcs(j)%active
    s_rdwbo = env%calc%calcs(j)%rdwbo
    s_rddip = env%calc%calcs(j)%rddip
    s_rdqat = env%calc%calcs(j)%rdqat
    s_rdgrad = env%calc%calcs(j)%rdgrad
    s_apiclean = env%calc%calcs(j)%apiclean
    if (allocated(env%calc%calcs(j)%calcspace)) s_calcspace = env%calc%calcs(j)%calcspace
    if (allocated(env%calc%calcs(j)%solvmodel)) s_solvmodel = env%calc%calcs(j)%solvmodel
    if (allocated(env%calc%calcs(j)%solvent)) s_solvent = env%calc%calcs(j)%solvent
    if (allocated(env%calc%calcs(j)%draco)) s_draco = env%calc%calcs(j)%draco
    !>--- also keep the OLD method definition so we can roll back if the new
    !>--- level string is not one the modern calculator understands
    s_id = env%calc%calcs(j)%id
    if (allocated(env%calc%calcs(j)%other)) s_other = env%calc%calcs(j)%other
    if (allocated(env%calc%calcs(j)%binary)) s_binary = env%calc%calcs(j)%binary

    !>--- re-create the level. NOTE: create() calls deallocate() first, which
    !>--- also drops the cached %calcfile/%gradfile/%systemcall -- exactly what
    !>--- we want, because those bake in the old binary and scratch paths.
    !>--- iostat is passed on purpose: create() is FATAL on an unknown level
    !>--- string unless the caller says it can recover, and this caller can.
    call env%calc%calcs(j)%create(trim(env%gfnver),iostat=io_lvl)

    !>--- ROLLBACK GUARD: a level string create() does not know (e.g. the legacy
    !>--- composite flags '--gfn2@gff') leaves id = jobtype%unknown = 0, i.e. a
    !>--- dead calculator. In that case restore the previous method rather than
    !>--- break the run.
    if (io_lvl /= 0 .or. env%calc%calcs(j)%id == jobtype%unknown) then
      env%calc%calcs(j)%id = s_id
      if (allocated(s_other)) env%calc%calcs(j)%other = s_other
      if (allocated(s_binary)) env%calc%calcs(j)%binary = s_binary
      call env%calc%calcs(j)%autocomplete(j)
      write (stdout,'(/,a,1x,a,1x,a)') '> WARNING: method',trim(env%gfnver), &
      & 'is not available in the new calculator routines;'
      write (stdout,'(a)') '>          keeping the previously selected level for this step.'
    end if

    !>--- restore
    env%calc%calcs(j)%chrg = s_chrg
    env%calc%calcs(j)%uhf = s_uhf
    env%calc%calcs(j)%refine_lvl = s_refine_lvl
    env%calc%calcs(j)%weight = s_weight
    env%calc%calcs(j)%active = s_active
    env%calc%calcs(j)%rdwbo = s_rdwbo
    env%calc%calcs(j)%rddip = s_rddip
    env%calc%calcs(j)%rdqat = s_rdqat
    env%calc%calcs(j)%rdgrad = s_rdgrad
    env%calc%calcs(j)%apiclean = s_apiclean
    if (allocated(s_calcspace)) env%calc%calcs(j)%calcspace = s_calcspace
    if (allocated(s_solvmodel)) env%calc%calcs(j)%solvmodel = s_solvmodel
    if (allocated(s_solvent)) env%calc%calcs(j)%solvent = s_solvent
    if (allocated(s_draco)) env%calc%calcs(j)%draco = s_draco
    !>--- the xtb subprocess route must keep pointing at the user's binary
    if (env%calc%calcs(j)%id == jobtype%xtbsys) then
      env%calc%calcs(j)%binary = trim(env%ProgName)
    end if

    if (allocated(s_calcspace)) deallocate (s_calcspace)
    if (allocated(s_solvmodel)) deallocate (s_solvmodel)
    if (allocated(s_solvent)) deallocate (s_solvent)
    if (allocated(s_draco)) deallocate (s_draco)
    if (allocated(s_other)) deallocate (s_other)
    if (allocated(s_binary)) deallocate (s_binary)
  end do

  return
end subroutine env_apply_calclevel

!================================================================================!
subroutine confscript2i(env,tim)
  use iso_fortran_env,only:wp => real64
  use crest_data
  implicit none
  type(systemdata) :: env
  type(timer)   :: tim
  if (env%legacy) then
    call confscript2i_legacy(env,tim)
  else
    if (.not.env%entropic) then
      call crest_search_imtdgc(env,tim)
    else
      call crest_search_entropy(env,tim)
    end if
  end if
end subroutine confscript2i

!================================================================================!
subroutine mdopt(env,tim)
  use iso_fortran_env,only:wp => real64
  use crest_data
  implicit none
  type(systemdata) :: env
  type(timer)   :: tim
  if (env%legacy) then
    call mdopt_legacy(env,tim)
  else
    call crest_ensemble_optimization(env,tim)
  end if
end subroutine mdopt

!================================================================================!
subroutine screen(env,tim)
  use iso_fortran_env,only:wp => real64
  use crest_data
  implicit none
  type(systemdata) :: env
  type(timer)   :: tim
  if (env%legacy) then
    call screen_legacy(env,tim)
  else
    call crest_ensemble_screening(env,tim)
  end if
end subroutine screen

!=================================================================================!

subroutine xtbsp(env,xtblevel)
  use iso_fortran_env,only:wp => real64
  use crest_data
  use strucrd,only:coord
  implicit none
  type(systemdata) :: env
  integer,intent(in),optional :: xtblevel
  interface
    subroutine crest_xtbsp(env,xtblevel,molin)
      import :: systemdata,coord
      type(systemdata) :: env
      integer,intent(in),optional :: xtblevel
      type(coord),intent(in),optional :: molin
    end subroutine crest_xtbsp
  end interface
  if (env%legacy) then
    call xtbsp_legacy(env,xtblevel)
  else
    call crest_xtbsp(env,xtblevel)
  end if
end subroutine xtbsp
subroutine xtbsp2(fname,env)
  use iso_fortran_env,only:wp => real64
  use crest_data
  use strucrd
  implicit none
  type(systemdata) :: env
  character(len=*),intent(in) :: fname
  type(coord) :: mol
  interface
    subroutine crest_xtbsp(env,xtblevel,molin)
      import :: systemdata,coord
      type(systemdata) :: env
      integer,intent(in),optional :: xtblevel
      type(coord),intent(in),optional :: molin
    end subroutine crest_xtbsp
  end interface
  if (env%legacy) then
    call xtbsp2_legacy(fname,env)
  else
    call mol%open(trim(fname))
    call crest_xtbsp(env,xtblevel=-1,molin=mol)
  end if
end subroutine xtbsp2

!=================================================================================!

subroutine confscript1(env,tim)
  use crest_parameters
  use crest_data
  implicit none
  type(systemdata) :: env
  type(timer)   :: tim
  write (stdout,*)
  write (stdout,*) 'This runtype has been entirely deprecated.'
  write (stdout,*) 'You may try an older version of the program if you want to use it.'
  stop
end subroutine confscript1

!=================================================================================!

subroutine nciflexi(env,flexval)
  use crest_parameters
  use crest_data
  use strucrd
  implicit none
  type(systemdata) :: env
  type(coord) ::mol
  real(wp) :: flexval
  if (env%legacy) then
    call nciflexi_legacy(env,flexval)
  else
    call env%ref%to(mol)
    call nciflexi_gfnff(mol,flexval)
  end if
end subroutine nciflexi

!================================================================================!

subroutine thermo_wrap(env,pr,nat,at,xyz,dirname, &
        &  nt,temps,et,ht,gt,stot,bhess)
!******************************************
!* Wrapper for a Hessian calculation
!* to get thermodynamics of the molecule
!*****************************************
  use crest_parameters,only:wp
  use crest_data
  implicit none
  !> INPUT
  type(systemdata) :: env
  logical,intent(in) :: pr
  integer,intent(in) :: nat
  integer,intent(inout) :: at(nat)
  real(wp),intent(inout) :: xyz(3,nat)  !> in Angstroem!
  character(len=*) :: dirname
  integer,intent(in)  :: nt
  real(wp),intent(in)  :: temps(nt)
  logical,intent(in) :: bhess       !> calculate bhess instead?
  !> OUTPUT
  real(wp),intent(out) :: et(nt)    !> enthalpy in Eh
  real(wp),intent(out) :: ht(nt)    !> enthalpy in Eh
  real(wp),intent(out) :: gt(nt)    !> free energy in Eh
  real(wp),intent(out) :: stot(nt)  !> entropy in cal/molK

  if (env%legacy) then
    call thermo_wrap_legacy(env,pr,nat,at,xyz,dirname, &
    &                    nt,temps,et,ht,gt,stot,bhess)
  else
    call thermo_wrap_new(env,pr,nat,at,xyz,dirname, &
    &                    nt,temps,et,ht,gt,stot,bhess)
  end if
end subroutine thermo_wrap

!================================================================================!

subroutine trialMD(env)
!***********************************************
!* subroutine trialMD
!* Takes the global metadynamics settings
!* And performs a short 1 ps simulation to
!* check if the molecular dynamis/metadynamics
!* will run, or if the timestep is too large
!***********************************************
  use crest_parameters,only:wp
  use crest_data
  implicit none
  !> INPUT
  type(systemdata) :: env

  if (env%legacy) then
    !> old xtb subprocess version
    call trialMD_legacy(env)
  else
    !> new calculator implementation
    call trialMD_calculator(env)
  end if

end subroutine trialMD

!================================================================================!

subroutine trialOPT(env)
!**********************************************************
!* subroutine trialOPT
!* Performs a geometry optimization of the structure
!* saved to env%ref and checks for changes in the topology
!**********************************************************
  use crest_data
  use crest_parameters,only:stdout
  implicit none
  !> INPUT
  type(systemdata) :: env

  if (env%legacy) then
    call xtbopt_legacy(env)
  else
    call trialOPT_calculator(env)
  end if

  if (env%crestver == crest_trialopt) then
!>-- if we reach this point in the standalone trialopt the geometry is ok!
    write (stdout,*)
    stop 'Geometry ok!'
  end if
end subroutine trialOPT

!================================================================================!

subroutine protonate(env,tim)
!*****************************************************
!* subroutine protonate
!* driver for the automated protonation site search
!* originally published in JCC
!*****************************************************
  use iso_fortran_env,only:wp => real64
  use crest_data
  implicit none
  type(systemdata) :: env
  type(timer)   :: tim
  if (env%legacy) then
    call protonate_legacy(env,tim)
  else
    call crest_new_protonate(env,tim)
  end if
end subroutine protonate

!================================================================================!

subroutine deprotonate(env,tim)
!*****************************************************
!* subroutine deprotonate
!* driver for the automated deprotonation site search
!*****************************************************
  use iso_fortran_env,only:wp => real64
  use crest_data
  implicit none
  type(systemdata) :: env
  type(timer)   :: tim
  if (env%legacy) then
    call deprotonate_legacy(env,tim)
  else
    call crest_new_deprotonate(env,tim)
  end if
end subroutine deprotonate

!================================================================================!

subroutine tautomerize(env,tim)
!*****************************************************
!* subroutine tautomerize
!* driver for the automated tautomer search
!*****************************************************
  use iso_fortran_env,only:wp => real64
  use crest_data
  implicit none
  type(systemdata) :: env
  type(timer)   :: tim
  if (env%legacy) then
    call tautomerize_legacy(env,tim)
  else
    call crest_new_tautomerize(env,tim)
  end if
end subroutine tautomerize

!========================================================================================!

subroutine catchdiatomic(env)
!****************************************
!* subroutine catchdiatomic
!* if we only have one or two atoms just
!* write the "optimized" structure
!****************************************
  use crest_data
  use crest_parameters
  use strucrd
  use iomod,only:copy
  use cregen_interface
  implicit none
  type(systemdata) :: env
  integer :: ich
  type(coord) :: mol
  if (env%legacy) then
    call catchdiatomic_legacy(env)
  else
    call env%ref%to(mol)
    open (file=conformerfile,newunit=ich)
    mol%xyz = mol%xyz*bohr !to ang
    call wrxyz(ich,mol%nat,mol%at,mol%xyz,env%ref%etot)
    close (ich)
    call copy('xtbopt.xyz',conformerfile)
    call copy(conformerfile,'crest_rotamers.xyz')
    call copy(conformerfile,'crest_best.xyz')
  end if
  call newcregen(env,6,'crest_rotamers.xyz')
end subroutine catchdiatomic
