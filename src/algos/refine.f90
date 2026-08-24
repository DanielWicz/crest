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

subroutine crest_refine(env,input,output)
!*******************************************************************
!* subroutine crest_refine
!* This subroutine will process the ensemble specified by "input"
!* and either overwrite it or write a new file "output".
!* The routine is intended to be called after geometry optimization
!* to re-rank structures by singlepoint energies, or add other
!* contributions to the energy.
!*******************************************************************
  use crest_parameters
  use crest_data
  use crest_calculator
  use strucrd
  use cregen_interface
  use parallel_interface
  implicit none
  type(systemdata),intent(inout) :: env
  character(len=*),intent(in) :: input
  character(len=*),intent(in),optional :: output
!===========================================================!
  integer :: i,j,k,l,io,ich,m,t1,t2
  logical :: pr,wr,ex
!===========================================================!
  character(len=:),allocatable :: outname
  real(wp) :: energy,gnorm
  real(wp),allocatable :: grad(:,:)
  character(len=:),allocatable :: ensnam
  integer :: nat,nall
  real(wp),allocatable :: eread(:),etmp(:)
  real(wp),allocatable :: xyz(:,:,:)
  integer,allocatable  :: at(:)
  integer :: nrefine,refine_stage
  integer :: nfail,jj
  logical,allocatable :: okmask(:)
  !> Energy handed to a structure whose refinement failed. This is the
  !> convention the rest of CREST already uses for failed optimizations
  !> (see src/algos/optimization.f90, "failed optimizations are assigned an
  !> energy of +1.0"): any positive value is hundreds of kcal/mol above every
  !> real total energy, so CREGEN's energy window removes it on the next sort.
  real(wp),parameter :: efail = 1.0_wp
!===========================================================!
!>--- setup
  if (present(output)) then
    outname = output !> new file
  else
    outname = input  !> overwrite
  end if
 
!>--- presorting step, if necessary
  if(env%refine_presort)then
    call newcregen(env,0,input)
    call rename('crest_ensemble.xyz',input)
  endif

!>--- read in
  call rdensemble(input,nat,nall,at,xyz,eread)
  allocate (etmp(nall),source=0.0_wp)
  allocate (okmask(nall),source=.true.)
!>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!>--- Important: crest_sploop requires coordinates in Bohrs
    xyz = xyz / bohr
!>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<!

!===========================================================!
  DO_REFINE: if (allocated(env%refine_queue)) then
!===========================================================!

    call smallhead('ensemble refinement')

    nrefine = size(env%refine_queue,1)

    do i = 1,nrefine
      refine_stage = env%refine_queue(i)
      !> set the calculator to the correct stage
      env%calc%refine_stage = refine_stage

      !>--- A refinement stage with no live calculator behind it would quietly
      !>--- re-rank the whole ensemble at whatever potential_core returns for an
      !>--- unimplemented job type. Refuse instead.
      if (refine_stage /= refine%confsolv) then
        if (.not.refine_stage_is_live(env,refine_stage)) then
          write (stderr,'(/,a,1x,i0,1x,a)') '**ERROR** refinement stage', &
          & refine_stage,'has no usable calculation level behind it.'
          write (stderr,'(a)') '          Refusing to re-rank the ensemble with it.'
          flush (stderr)
          env%iostatus_meta = status_failed
          return
        end if
      end if

      select case (refine_stage)
      case (refine%singlepoint)
        write (stdout,'("> Singlepoint re-ranking for ",i0," structures")') nall
        call crest_sploop(env,nat,nall,at,xyz,eread,ok=okmask)
        !> A failed singlepoint leaves eread = 0.0 Eh, which is NOT "no energy":
        !> it is above every real total energy and would look like a legitimate,
        !> very high-lying structure. Mark it with the standard failure energy.
        call mark_failed(nall,eread,okmask,efail,nfail)
        if (nfail > 0) write (stdout,'("> ",i0,a)') nfail, &
        & ' structure(s) failed and were flagged for removal'

      case (refine%correction)
        write (stdout,'("> Additive correction for ",i0," structures")') nall
        call crest_sploop(env,nat,nall,at,xyz,etmp,ok=okmask)
        !> PHYSICS: this branch is the dangerous one. It used to add etmp
        !> unconditionally, and a failed structure contributed etmp = 0.0, so it
        !> stayed on the UNCORRECTED level while every competitor moved to the
        !> corrected one -- two different Hamiltonians ranked against each other
        !> in a single list. Flag such structures instead of mixing levels.
        do jj = 1,nall
          if (okmask(jj)) then
            eread(jj) = eread(jj)+etmp(jj)
          else
            eread(jj) = efail
          end if
        end do
        nfail = count(.not.okmask)
        if (nfail > 0) write (stdout,'("> ",i0,a)') nfail, &
        & ' structure(s) could not be corrected and were flagged for removal'

      case (refine%geoopt)
        write (stdout,'("> Geometry optimization of ",i0," structures")') nall
        !> crest_oloop already assigns its own +1.0 failure energy and leaves
        !> the input geometry in place for structures that did not optimize.
        call crest_oloop(env,nat,nall,at,xyz,eread,.false.)

      case(refine%confsolv)
        call new_ompautoset(env,'subprocess',1,t1,t2)
        write (stdout,'("> ConfSolv: ΔΔGsoln estimation from 3D directed message passing neural networks (D-MPNN)")')
        call confsolv_request( input, nall, t2, etmp, io)
        if(io == 0)then
        eread(:) = etmp(:)*kcaltoau  !> since CREGEN deals with Eh energies
        else
        !> Silently keeping the old energies here would report a solution-phase
        !> ranking that is really the gas-phase one. Fail instead.
        write (stderr,'(/,a)') '**ERROR** ConfSolv request failed; the ensemble was NOT rescored.'
        flush (stderr)
        env%iostatus_meta = status_failed
        return
        endif
      end select
      write(stdout,*) 
    end do

    !> reset the refinement stage of the calculator
    env%calc%refine_stage = 0

!===========================================================!
  end if DO_REFINE
!===========================================================!

!>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!>--- Important: ensemble file must be written in AA
  xyz = xyz / angstrom
!>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!>--- write output ensemble
  call wrensemble(outname,nat,nall,at,xyz,eread)

!===========================================================!
  deallocate (etmp,eread,xyz,at)
  if (allocated(okmask)) deallocate (okmask)
  return

contains

  subroutine mark_failed(n,e,ok,ef,nf)
!*******************************************************
!* Replace the energies of failed structures with the
!* standard failure energy so CREGEN drops them, instead
!* of leaving a value that looks like a real energy.
!*******************************************************
    implicit none
    integer,intent(in) :: n
    real(wp),intent(inout) :: e(n)
    logical,intent(in) :: ok(n)
    real(wp),intent(in) :: ef
    integer,intent(out) :: nf
    integer :: m
    nf = 0
    do m = 1,n
      if (.not.ok(m)) then
        e(m) = ef
        nf = nf+1
      end if
    end do
  end subroutine mark_failed

  logical function refine_stage_is_live(env,stage) result(live)
!*******************************************************
!* Is there an active calculation level tagged with this
!* refinement stage, and does it have a real job type?
!*******************************************************
    implicit none
    type(systemdata),intent(in) :: env
    integer,intent(in) :: stage
    integer :: m
    live = .false.
    do m = 1,env%calc%ncalculations
      if (env%calc%calcs(m)%refine_lvl /= stage) cycle
      if (.not.env%calc%calcs(m)%active) cycle
      if (env%calc%calcs(m)%id == jobtype%unknown) cycle
      !> A geometry optimization needs forces. A level configured not to read a
      !> gradient would "converge" at the input geometry with a high-level
      !> energy, which looks like a successful re-optimization and is not one.
      if (stage == refine%geoopt .and. .not.env%calc%calcs(m)%rdgrad) cycle
      live = .true.
      return
    end do
  end function refine_stage_is_live

end subroutine crest_refine
!========================================================================================!
!========================================================================================!
