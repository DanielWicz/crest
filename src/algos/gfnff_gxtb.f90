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
  use crest_parameters,only:wp,stdout,stderr,autokcal,bohr
  use crest_data
  use strucrd
  use iomod
  use utilities
  use cregen_interface
  use gxtb_reopt_module
  use miscdata,only:rcov
  use crest_cn_module
  implicit none
  private

  public :: crest_gfnff_gxtb_stage
  public :: crest_gxtb_replay
  public :: gfnff_gxtb_stage_count
  public :: gxtb_flag_failure
  public :: gxtb_clear_failure_marker
  public :: gxtb_write_populations

  !> the stored input of the post-sampling stage (GFN-FF energies)
  character(len=*),parameter :: stored_ensemble = 'crest_gfnff_ensemble.xyz'

  !> how often the stage was entered in this program run (queried by crest_main,
  !> which insists on exactly one run for every -gfnff:gxtb workflow)
  integer,save :: nstage_runs = 0
  character(len=*),parameter :: failure_marker = 'CREST_GXTB_FAILED'

  !> set when THIS run is about to write (or clear) crest_conformers.xyz and the
  !> other result files; only then may a failure remove them (never in -nocregen)
  logical,save :: run_writes_results = .false.

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
    integer :: cs,gs
    character(len=256) :: msg

    nstage_runs = nstage_runs+1
    !>--- 1. keep the GFN-FF ensemble as the stored input of this stage
    call execute_command_line("cp -f '"//ensemble_file//"' '"//stored_ensemble//"'",exitstat=io,cmdstat=cs,cmdmsg=msg)
    if (io /= 0 .or. cs /= 0) then
      write (stderr,'(a)') '**ERROR** could not store '//stored_ensemble
      env%iostatus_meta = status_failed
      return
    end if
    call rdensembleparam(ensemble_file,nat,nin)

    !>--- 2. presort: CREGEN at GFN-FF energies, unique conformers only
    run_writes_results = .true.
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
    call getcwd(cdir,gs)
    cwd = trim(cdir)
    if (gs /= 0) then
      write (stderr,'(a)') '**ERROR** cannot determine the working directory.'
      env%iostatus_meta = status_failed
      return
    end if
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
    call reopt%cleanup()

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
    integer :: nat,nall,io,cs
    character(len=256) :: msg
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
    if (env%gxtb_nocregen) then
      call reoptimize_given_structures(env,nat,nall)
      call tim%stop(14)
      return
    end if
    write (stdout,*)
    write (stdout,'(3x,a)') '====================================================='
    write (stdout,'(3x,a)') '|  g-xTB gas-phase re-optimization (replay)          |'
    write (stdout,'(3x,a)') '====================================================='
    write (stdout,'(1x,a,a)') 'Ensemble: ',trim(env%ensemblename)

    !>--- clean start: previous results must not influence this run
    run_writes_results = .true.
    call rmcres()
    call checkname_xyz(crefile,inpnam,outnam)
    call execute_command_line("cp -f '"//trim(env%ensemblename)//"' '"//trim(outnam)//"'",exitstat=io,cmdstat=cs,cmdmsg=msg)
    if (io /= 0 .or. cs /= 0) then
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

!================================================================================!

  subroutine reoptimize_given_structures(env,nat_in,nall_in)
    !> "--gxtbreopt <file> -nocregen": no CREGEN of any kind, no energy window,
    !> no merging, no topology filter. Every structure of the file, in file
    !> order, is re-optimized; crest_gxtb_reopt.xyz holds one block per INPUT
    !> index. This path counts as the one run of the stage (nstage_runs), so
    !> the "ran exactly once" assertion of crest_main holds here as well.
    type(systemdata),intent(inout) :: env
    integer,intent(in) :: nat_in,nall_in
    integer :: nat,nall
    character(len=*),parameter :: outfile = 'crest_gxtb_reopt.xyz'
    type(gxtb_reoptimizer) :: reopt
    integer,allocatable :: at(:),idx(:)
    real(wp),allocatable :: xyz(:,:,:),eread(:)
    logical,allocatable :: ok(:)
    integer,allocatable :: topo(:)
    logical :: good,written
    integer :: t0,t1,rate,k,gs
    character(len=512) :: cdir
    character(len=:),allocatable :: cwd,line
    character(len=32) :: ctmp

    nstage_runs = nstage_runs+1
    nat = nat_in
    nall = nall_in
    call rdensemble(env%ensemblename,nat,nall,at,xyz,eread)
    allocate (ok(nall))
    call getcwd(cdir,gs)
    cwd = trim(cdir)
    if (gs /= 0) then
      write (stderr,'(a)') '**ERROR** cannot determine the working directory.'
      env%iostatus_meta = status_failed
      return
    end if
    call reopt%init(trim(env%ProgName),env%chrg,env%uhf,env%gxtb_optlevel,max(1,env%threads), &
    &    cwd//'/.gxtb_reopt',env%keepModef,good)
    if (.not.good) then
      write (stderr,'(a)') '**ERROR** invalid g-xTB re-optimization settings.'
      env%iostatus_meta = status_failed
      return
    end if
    write (stdout,'(/,1x,a)') 'g-xTB gas-phase re-optimization of the given structures (no CREGEN)'
    call system_clock(t0,rate)
    call reopt%reoptimize(nat,nall,at,xyz,eread,ok)
    call system_clock(t1)
    allocate (topo(nall))
    call connectivity_flags(env,nat,nall,at,xyz,ok,topo)
    call reopt%write_results(outfile,nat,nall,ok,topo,written)
    call reopt%cleanup()
    idx = reopt%failed_indices()
    write (stdout,'(3x,a,i0)') 'structures in the input file      : ',nall
    write (stdout,'(3x,a,i0)') 'successfully re-optimized         : ',count(ok)
    write (stdout,'(3x,a,i0)') 'failed (ok=0 in the output)       : ',reopt%nfailed()
    if (size(idx) > 0) then
      line = ''
      do k = 1,size(idx)
        write (ctmp,'(i0)') idx(k)
        line = line//' '//trim(ctmp)
      end do
      write (stdout,'(3x,a)') 'failed structure indices (1-based)  :'//line
    end if
    write (stdout,'(3x,a,f10.2,a)') 'wall time                         : ', &
    &    real(t1-t0,wp)/real(rate,wp),' s'
    write (stdout,'(3x,a)') 'results written to '//outfile
    if (.not.written) then
      write (stderr,'(a)') '**ERROR** could not write '//outfile
      env%iostatus_meta = status_failed
    end if
  end subroutine reoptimize_given_structures

!================================================================================!

  subroutine connectivity_flags(env,nat,nall,at,xyz_ang,ok,topo)
    !> topo(k) = 1 if structure k has the same connectivity as the reference
    !> structure ("coord", i.e. the input structure) by the criterion of
    !> CREGEN's topology check (cregen_topocheck: calc_ncoord + bondtotopo,
    !> same neighbour-list comparison), 0 if not. Failed structures get -1.
    type(systemdata),intent(in) :: env
    integer,intent(in) :: nat,nall
    integer,intent(in) :: at(nat)
    real(wp),intent(in) :: xyz_ang(3,nat,nall)
    logical,intent(in) :: ok(nall)
    integer,intent(out) :: topo(nall)
    integer :: j,ntopo
    integer,allocatable :: atdum(:),toporef(:),topoj(:)
    real(wp),allocatable :: cref(:,:),c1(:,:),cn(:),bond(:,:)
    logical,allocatable :: neighmat(:,:)
    topo = -1
    ntopo = nat*(nat+1)/2
    allocate (cref(3,nat),atdum(nat),c1(3,nat),cn(nat),bond(nat,nat))
    allocate (toporef(ntopo),topoj(ntopo),neighmat(nat,nat))
    neighmat = .false.
    call rdcoord('coord',nat,atdum,cref)
    cn = 0.0_wp
    bond = 0.0_wp
    call calc_ncoord(nat,atdum,cref,rcov,cn,400.0_wp,bond)
    if (allocated(env%excludeTOPO)) then
      call bondtotopo_excl(nat,at,bond,cn,ntopo,toporef,neighmat,env%excludeTOPO)
    else
      call bondtotopo(nat,at,bond,cn,ntopo,toporef,neighmat)
    end if
    do j = 1,nall
      if (.not.ok(j)) cycle
      c1(:,:) = xyz_ang(:,:,j)/bohr
      cn = 0.0_wp
      bond = 0.0_wp
      call calc_ncoord(nat,at,c1,rcov,cn,400.0_wp,bond)
      if (allocated(env%excludeTOPO)) then
        call bondtotopo_excl(nat,at,bond,cn,ntopo,topoj,neighmat,env%excludeTOPO)
      else
        call bondtotopo(nat,at,bond,cn,ntopo,topoj,neighmat)
      end if
      topo(j) = merge(1,0,all(toporef == topoj))
    end do
  end subroutine connectivity_flags

!================================================================================!

  subroutine gxtb_flag_failure(reason)
    !> A -gfnff:gxtb workflow failed: remove the GFN-FF-level result files that
    !> could be mistaken for final results (exact names) and leave a marker.
    character(len=*),intent(in) :: reason
    integer :: ich
    !> only files of this run: a caller's crest_conformers.xyz etc. must survive
    if (run_writes_results) then
      call rmcres()   !> crest_conformers.xyz and crest_rotamers*.xyz (this run's working copies)
      call remove('crest_best.xyz')
      call remove('crest.energies')
    end if
    open (newunit=ich,file=failure_marker,status='replace')
    write (ich,'(a)') 'The -gfnff:gxtb / --gxtbreopt workflow FAILED: '//trim(reason)
    write (ich,'(a)') 'Result files written by this run at the GFN-FF level were removed; this run produced no g-xTB result.'
    close (ich)
    write (stderr,'(a)') '**ERROR** gxtb workflow failed ('//trim(reason)//'); see '//failure_marker
  end subroutine gxtb_flag_failure

  subroutine gxtb_clear_failure_marker()
    call remove(failure_marker)
  end subroutine gxtb_clear_failure_marker

!================================================================================!

  subroutine gxtb_write_populations(env)
    !> Plain Boltzmann populations of the conformers in crest_conformers.xyz
    !> (no rotamer degeneracy), to crest_gxtb_populations.dat and stdout.
    type(systemdata),intent(in) :: env
    character(len=*),parameter :: pfile = 'crest_gxtb_populations.dat'
    real(wp),parameter :: kb = 3.166811563e-6_wp   !> Eh/K
    integer :: nat,nall,i,ich,io
    integer,allocatable :: at(:)
    real(wp),allocatable :: xyz(:,:,:),e(:),w(:)
    real(wp) :: kt,s
    logical :: ex
    inquire (file=conformerfile,exist=ex)
    if (.not.ex) return
    call rdensembleparam(conformerfile,nat,nall)
    if (nall < 1) return
    call rdensemble(conformerfile,nat,nall,at,xyz,e)
    kt = kb*env%tboltz
    allocate (w(nall))
    w = exp(-(e-e(1))/kt)
    s = sum(w)
    w = w/s
    open (newunit=ich,file=pfile,status='replace',iostat=io)
    if (io /= 0) return
    call table(ich)
    close (ich)
    call table(stdout)
  contains
    subroutine table(ch)
      integer,intent(in) :: ch
      write (ch,'(a)') '# g-xTB conformer populations: plain Boltzmann weights, NO rotamer degeneracy.'
      write (ch,'(a,es22.15,a,f8.3,a)') '# kB = ',kb,' Eh/K, T = ',env%tboltz,' K; Erel in kcal/mol'
      write (ch,'(a)') '# (CREGEN populations sum rotamer weights and use a rounded kT; they are not these.)'
      write (ch,'(a)') '# index   E(Eh)               Erel(kcal/mol)    weight'
      do i = 1,nall
        write (ch,'(i6,2x,f18.8,2x,f12.6,2x,es24.16)') i,e(i),(e(i)-e(1))*autokcal,w(i)
      end do
    end subroutine table
  end subroutine gxtb_write_populations

end module gfnff_gxtb_module
