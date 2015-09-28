!
! Copyright (C) 2003-2013 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE vloc_psi_gamma(lda, n, m, psi, v, hpsi)
  !-----------------------------------------------------------------------
  !
  ! Calculation of Vloc*psi using dual-space technique - Gamma point
  !
  USE parallel_include
  USE kinds,   ONLY : DP
  USE gvecs, ONLY : nls, nlsm
  USE wvfct,   ONLY : igk
  USE mp_bands,      ONLY : me_bgrp
  USE fft_base,      ONLY : dffts, tg_gather
  USE fft_interfaces,ONLY : fwfft, invfft
  USE wavefunctions_module,  ONLY: psic
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: lda, n, m
  COMPLEX(DP), INTENT(in)   :: psi (lda, m)
  COMPLEX(DP), INTENT(inout):: hpsi (lda, m)
  REAL(DP), INTENT(in) :: v(dffts%nnr)
  !
  INTEGER :: ibnd, j, incr
  COMPLEX(DP) :: fp, fm
  !
  LOGICAL :: use_tg
  ! Variables for task groups
  REAL(DP),    ALLOCATABLE :: tg_v(:)
  COMPLEX(DP), ALLOCATABLE :: tg_psic(:)
  INTEGER :: v_siz, idx, ioff
  !
#if defined(__CUDA) && !defined(__DISABLE_CUDA_VLOCPSI) && ( !defined(__PARA) || defined(__USE_3D_FFT) )
  CALL vloc_psi_gamma_gpu ( lda, n, m, psi, v, hpsi )
  RETURN
#endif
  !
  incr = 2
  !
  ! The following is dirty trick to prevent usage of task groups if
  ! the number of bands is smaller than the number of task groups 
  ! 
  use_tg = dffts%have_task_groups
  dffts%have_task_groups  = dffts%have_task_groups .and. ( m >= dffts%nogrp )
  !
  IF( dffts%have_task_groups ) THEN
     !
     v_siz =  dffts%tg_nnr * dffts%nogrp
     !
     ALLOCATE( tg_v   ( v_siz ) )
     ALLOCATE( tg_psic( v_siz ) )
     !
     CALL tg_gather( dffts, v, tg_v )
     !
     incr = 2 * dffts%nogrp
     !
  ENDIF
  !
  ! the local potential V_Loc psi. First bring psi to real space
  !
  DO ibnd = 1, m, incr
     !
     IF( dffts%have_task_groups ) THEN
        !
        tg_psic = (0.d0, 0.d0)
        ioff   = 0
        !
        DO idx = 1, 2*dffts%nogrp, 2
           IF( idx + ibnd - 1 < m ) THEN
              DO j = 1, n
                 tg_psic(nls (igk(j))+ioff) =        psi(j,idx+ibnd-1) + &
                                      (0.0d0,1.d0) * psi(j,idx+ibnd)
                 tg_psic(nlsm(igk(j))+ioff) = conjg( psi(j,idx+ibnd-1) - &
                                      (0.0d0,1.d0) * psi(j,idx+ibnd) )
              ENDDO
           ELSEIF( idx + ibnd - 1 == m ) THEN
              DO j = 1, n
                 tg_psic(nls (igk(j))+ioff) =        psi(j,idx+ibnd-1)
                 tg_psic(nlsm(igk(j))+ioff) = conjg( psi(j,idx+ibnd-1) )
              ENDDO
           ENDIF

           ioff = ioff + dffts%tg_nnr

        ENDDO
        !
     ELSE
        !
        psic(:) = (0.d0, 0.d0)
        IF (ibnd < m) THEN
           ! two ffts at the same time
           DO j = 1, n
              psic(nls (igk(j)))=      psi(j,ibnd) + (0.0d0,1.d0)*psi(j,ibnd+1)
              psic(nlsm(igk(j)))=conjg(psi(j,ibnd) - (0.0d0,1.d0)*psi(j,ibnd+1))
           ENDDO
        ELSE
           DO j = 1, n
              psic (nls (igk(j))) =       psi(j, ibnd)
              psic (nlsm(igk(j))) = conjg(psi(j, ibnd))
           ENDDO
        ENDIF
        !
     ENDIF
     !
     !   fft to real space
     !   product with the potential v on the smooth grid
     !   back to reciprocal space
     !
     IF( dffts%have_task_groups ) THEN
        !
        CALL invfft ('Wave', tg_psic, dffts)
        !
        DO j = 1, dffts%nr1x*dffts%nr2x*dffts%tg_npp( me_bgrp + 1 )
           tg_psic (j) = tg_psic (j) * tg_v(j)
        ENDDO
        !
        CALL fwfft ('Wave', tg_psic, dffts)
        !
     ELSE
        !
        CALL invfft ('Wave', psic, dffts)
        !
        DO j = 1, dffts%nnr
           psic (j) = psic (j) * v(j)
        ENDDO
        !
        CALL fwfft ('Wave', psic, dffts)
        !
     ENDIF
     !
     !   addition to the total product
     !
     IF( dffts%have_task_groups ) THEN
        !
        ioff   = 0
        !
        DO idx = 1, 2*dffts%nogrp, 2
           !
           IF( idx + ibnd - 1 < m ) THEN
              DO j = 1, n
                 fp= ( tg_psic( nls(igk(j)) + ioff ) +  &
                       tg_psic( nlsm(igk(j)) + ioff ) ) * 0.5d0
                 fm= ( tg_psic( nls(igk(j)) + ioff ) -  &
                       tg_psic( nlsm(igk(j)) + ioff ) ) * 0.5d0
                 hpsi (j, ibnd+idx-1) = hpsi (j, ibnd+idx-1) + &
                                        cmplx( dble(fp), aimag(fm),kind=DP)
                 hpsi (j, ibnd+idx  ) = hpsi (j, ibnd+idx  ) + &
                                        cmplx(aimag(fp),- dble(fm),kind=DP)
              ENDDO
           ELSEIF( idx + ibnd - 1 == m ) THEN
              DO j = 1, n
                 hpsi (j, ibnd+idx-1) = hpsi (j, ibnd+idx-1) + &
                                         tg_psic( nls(igk(j)) + ioff )
              ENDDO
           ENDIF
           !
           ioff = ioff + dffts%nr3x * dffts%nsw( me_bgrp + 1 )
           !
        ENDDO
        !
     ELSE
        IF (ibnd < m) THEN
           ! two ffts at the same time
           DO j = 1, n
              fp = (psic (nls(igk(j))) + psic (nlsm(igk(j))))*0.5d0
              fm = (psic (nls(igk(j))) - psic (nlsm(igk(j))))*0.5d0
              hpsi (j, ibnd)   = hpsi (j, ibnd)   + &
                                 cmplx( dble(fp), aimag(fm),kind=DP)
              hpsi (j, ibnd+1) = hpsi (j, ibnd+1) + &
                                 cmplx(aimag(fp),- dble(fm),kind=DP)
           ENDDO
        ELSE
           DO j = 1, n
              hpsi (j, ibnd)   = hpsi (j, ibnd)   + psic (nls(igk(j)))
           ENDDO
        ENDIF
     ENDIF
     !
  ENDDO
  !
  IF( dffts%have_task_groups ) THEN
     !
     DEALLOCATE( tg_psic )
     DEALLOCATE( tg_v )
     !
  ENDIF
  dffts%have_task_groups = use_tg
  !
  RETURN
END SUBROUTINE vloc_psi_gamma
!

!dir$ attributes offload:mic :: mic_fftw_pseudopotentialmkl
subroutine mic_fftw_pseudopotentialmkl(nppm, ptr_cmp,v, nr1,nr2,nr1x,nr2x,nnr)
  USE kinds,   ONLY : DP
  USE ifport,  ONLY: dclock
  USE OMP_LIB
  use MKL_DFTI
  use iso_c_binding
  implicit none 
  double precision :: tscale
  integer :: nppm
  integer :: inca
  integer :: nr1,nr2,nr1x,nr2x,nnr
  COMPLEX(DP), dimension(nnr)  :: ptr_cmp
  REAL(DP) ,dimension(nnr) :: v
  integer :: j
  integer :: idir
  logical,save :: initplan2 = .true.
  integer :: status = 0 
  type(DFTI_DESCRIPTOR), POINTER, SAVE :: hand
  
  inca=1
  IF( initplan2 .eq. .true. ) THEN
     hand => null()
     status = DftiCreateDescriptor(hand, DFTI_DOUBLE, DFTI_COMPLEX, 2,(/nr1,nr2/))
     if(status /= 0) then 
        write(*,*) "stopped in DftiCreateDescriptor", status
        stop
     endif
     status = DftiSetValue(hand, DFTI_NUMBER_OF_TRANSFORMS,nppm)
     if(status /= 0)then 
        write(*,*) "stopped in DFTI_NUMBER_OF_TRANSFORMS", status
        stop
     endif
     status = DftiSetValue(hand,DFTI_INPUT_DISTANCE, nr1x*nr2x )
     if(status /= 0)then 
        write(*,*) "stopped in DFTI_INPUT_DISTANCE", status
        stop
     endif
     status = DftiSetValue(hand, DFTI_PLACEMENT, DFTI_INPLACE)     
     if(status /= 0)then 
        write(*,*) "stopped in DFTI_PLACEMENT", status
        stop
     endif
     tscale = 1.0_DP/ (nr1 * nr2 )
     Status = DftiSetValue( hand, DFTI_FORWARD_SCALE, tscale);
     if(status /= 0)then 
        write(*,*) "stopped in DFTI_FORWARD_SCALE", status
        stop
     endif
     Status = DftiSetValue( hand, DFTI_BACKWARD_SCALE, dble(1) );
     if(status /= 0)then 
        write(*,*) "stopped in DFTI_BACKWARD_SCALE", status
        stop
     endif
     status = DftiCommitDescriptor(hand)
     if(status /= 0)then 
        write(*,*) "stopped in DftiCommitDescriptor", status
        stop
     endif
     initplan2 = .false.
  END IF
  status = DftiComputeBackward(hand, ptr_cmp)
  if(status /= 0) then
     write(*,*) "stopped in DftiComputeBackward", status
     stop
  endif
  !
  !$omp parallel do
  DO j = 1, nnr
     ptr_cmp (j) = ptr_cmp (j) * v(j)
  ENDDO
  !$omp end parallel do
  !
  status = DftiComputeForward(hand, ptr_cmp)
  if(status /= 0) then 
     write(*,*) "stopped in DftiComputeForward", status
     stop
  endif
end subroutine mic_fftw_pseudopotentialmkl

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! THIS IS THE WORKING OFFLOAD VERSION
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!SUBROUTINE vloc_psi_k_catchwait(lda, n, m, psi, v, hpsi)
SUBROUTINE vloc_psi_k(lda, n, m, psi, v, hpsi)
  !-----------------------------------------------------------------------
  !
  ! Calculation of Vloc*psi using dual-space technique - k-points
  !
  USE parallel_include
  USE kinds,   ONLY : DP
  USE gvecs, ONLY : nls, nlsm
  USE wvfct,   ONLY : igk
  USE mp_bands,      ONLY : me_bgrp
  USE fft_base,      ONLY : dffts, tg_gather 
!  USE fft_interfaces,ONLY : fwfft_z, fwfft_xy,fwfft_scatter, invfft_z,invfft_xy,invfft_scatter
  USE fft_parallel,ONLY : tg_cft3s_scatter !, mic_fftw_pseudopotential
  USE ifport,  ONLY: dclock
  use MKL_DFTI
  USE wavefunctions_module,  ONLY: psic
  use, intrinsic :: iso_c_binding
  !
  IMPLICIT NONE
  !
!  include 'fftw3.f'
  INTEGER, INTENT(in) :: lda, n, m
  COMPLEX(DP), INTENT(in)   :: psi (lda, m)
  COMPLEX(DP), INTENT(inout):: hpsi (lda, m)
  REAL(DP), INTENT(in) :: v(dffts%nnr)
  !
  INTEGER :: ibnd, j, incr
  !
  LOGICAL :: use_tg
  ! Task Groups
  REAL(DP),    ALLOCATABLE :: tg_v(:)
  COMPLEX(DP), ALLOCATABLE :: tg_psic(:)
  INTEGER :: v_siz, idx, ioff
  !
  COMPLEX(DP), target, ALLOCATABLE :: psic_aux_1(:)
  COMPLEX(DP), target, ALLOCATABLE :: aux_1(:)
  COMPLEX(DP), target, ALLOCATABLE :: psic_aux_2(:)
  COMPLEX(DP), target, ALLOCATABLE :: aux_2(:)
  COMPLEX(DP), target, ALLOCATABLE :: psic_aux_3(:)
  COMPLEX(DP), target, ALLOCATABLE :: aux_3(:)
  COMPLEX(DP), pointer, dimension(:) :: ptr_snd
  COMPLEX(DP), pointer, dimension(:) :: ptr_snda
  COMPLEX(DP), pointer, dimension(:) :: ptr_cmp
  COMPLEX(DP), pointer, dimension(:) :: ptr_cmpa
  COMPLEX(DP), pointer, dimension(:) :: ptr_rec
  COMPLEX(DP), pointer, dimension(:) :: ptr_reca
  COMPLEX(DP), pointer, dimension(:) :: ptr_tmp
  LOGICAL, SAVE :: initplan = .true.
  integer :: idir, inca, nr1,nr2,nr1x,nr2x,nppm,nnr,nr3,nr3x,nswm
  integer*8 :: signal=7777
  integer :: status
  logical :: catchwait=.false.
  type(DFTI_DESCRIPTOR), POINTER, SAVE :: hand
  REAL(DP) :: tscale
  double precision :: t1 
  double precision :: t2
  character(len=32) :: carg
  integer, save :: offload_device=0

  t2=0
  t1=dclock()
  !
  !dir$ attributes offload:mic                :: psic_aux_1
  !dir$ attributes align: 4096                :: psic_aux_1
  !dir$ attributes offload:mic                :: aux_1
  !dir$ attributes align: 4096                :: aux_1
  !dir$ attributes offload:mic                :: psic_aux_2
  !dir$ attributes align: 4096                :: psic_aux_2
  !dir$ attributes offload:mic                :: aux_2
  !dir$ attributes align: 4096                :: aux_2
  !dir$ attributes offload:mic                :: psic_aux_3
  !dir$ attributes align: 4096                :: psic_aux_3
  !dir$ attributes offload:mic                :: aux_3
  !dir$ attributes align: 4096                :: aux_3
  !dir$ attributes offload:mic                :: mic_fftw_pseudopotentialmkl
  !
  ! -- for most purposes, we need parameters in a nice, bit-copyable format
  !
  IF( initplan .eq. .true. ) THEN
     call getenv("QE_MIC_DEVICE",carg)
     if(len(trim(carg)).ne.0) then
        read(carg,*) offload_device
     end if
  end IF

  inca = 1
  nppm = dffts%npp(dffts%mype+1)
  nr1 =  dffts%nr1
  nr2 =  dffts%nr2
  nr3 =  dffts%nr3
  nr1x = dffts%nr1x
  nr2x = dffts%nr2x 
  nr3x = dffts%nr3x 
  nnr = dffts%nnr
  nswm=dffts%nsw( dffts%mype + 1 )
!  write(*,*)"FFT n =", nr1,nr2,nr3
!  write(*,*)"FFT bands =", m,nnr
!  write(*,*)"FFT dim =", nppm, nr1x,nr2x

  !
  ! -- allocate the 3 fields and 3 auxiliary fields for the FFT runs - allocate the psic fields on the card
  allocate(psic_aux_1(dffts%nnr))
  !dir$ offload_transfer target(mic:offload_device) nocopy(psic_aux_1:length(nnr) alloc_if(.true.) free_if(.false.))
  allocate(aux_1(dffts%nnr))
  allocate(psic_aux_2(dffts%nnr))
  !dir$ offload_transfer target(mic:offload_device) nocopy(psic_aux_2:length(nnr) alloc_if(.true.) free_if(.false.))
  allocate(aux_2(dffts%nnr))
  allocate(psic_aux_3(dffts%nnr))
  !dir$ offload_transfer target(mic:offload_device) nocopy(psic_aux_3:length(nnr) alloc_if(.true.) free_if(.false.))
  allocate(aux_3(dffts%nnr))
  !
  ! -- unsure what happens if we have taskgroups - let's exclude that for now
  if(dffts%have_task_groups) then
     write(*,*) "No taskgroups, please!"
  end if
  !
  ! - establish an (arbitrary) relationship between pointers and buffers
  ptr_snd=>psic_aux_1
  ptr_snda=>aux_1
  ptr_cmp=>psic_aux_2
  ptr_cmpa=>aux_2
  ptr_rec=>psic_aux_3
  ptr_reca=>aux_3
  !
  !DIR$ offload_transfer target(mic:offload_device) in(v:length(nnr) alloc_if(.true.) free_if(.false.)) 
  !
  IF( initplan .eq. .true. ) THEN
     hand => null()
     status = DftiCreateDescriptor(hand, DFTI_DOUBLE, DFTI_COMPLEX, 1,nr3)
     if(status /= 0) then 
        write(*,*) "stopped in DftiCreateDescriptor", status
        stop
     endif
     status = DftiSetValue(hand, DFTI_NUMBER_OF_TRANSFORMS,nswm)
     if(status /= 0)then 
        write(*,*) "stopped in DFTI_NUMBER_OF_TRANSFORMS", status
        stop
     endif
     status = DftiSetValue(hand,DFTI_INPUT_DISTANCE, nr3x )
     if(status /= 0)then 
        write(*,*) "stopped in DFTI_INPUT_DISTANCE", status
        stop
     endif
     status = DftiSetValue(hand, DFTI_PLACEMENT, DFTI_INPLACE)     
     if(status /= 0)then 
        write(*,*) "stopped in DFTI_PLACEMENT", status
        stop
     endif
     tscale = 1.0_DP/nr3
     Status = DftiSetValue( hand, DFTI_FORWARD_SCALE, tscale);
     if(status /= 0)then 
        write(*,*) "stopped in DFTI_FORWARD_SCALE", status
        stop
     endif
     Status = DftiSetValue( hand, DFTI_BACKWARD_SCALE, dble(1) );
     if(status /= 0)then 
        write(*,*) "stopped in DFTI_BACKWARD_SCALE", status
        stop
     endif
     status = DftiCommitDescriptor(hand)
     if(status /= 0)then 
        write(*,*) "stopped in DftiCommitDescriptor", status
        stop
     endif
     
     initplan = .false.
  END IF
  ! -- make sure catchwait is false until set true
  catchwait=.false.
  DO ibnd = 0, m+1,1
     !
     if((ibnd.gt.0).and.(ibnd.lt.(m+1))) then
        catchwait=.true.
        ! -- this is the actual offloaded 2D FFT and pseudopotential application
        !DIR$ offload target(mic:offload_device) in(nppm, nr1, nr2, nr1x, nr2x, nnr) &
        !DIR$ & in(ptr_cmp:length(0) alloc_if(.false.) free_if(.false.)) &
        !DIR$ & in(v:length(0) alloc_if(.false.) free_if(.false.)) signal(signal)
        call  mic_fftw_pseudopotentialmkl(nppm, ptr_cmp,v, nr1,nr2,nr1x,nr2x,nnr)
        catchwait=.true.
     end if
     !
     ! -- begin the main loop 
     if(ibnd.lt.m) then
        !
        ptr_snd(:) = (0.d0, 0.d0)
        ptr_snd(nls (igk(1:n))) = psi(1:n,ibnd+1)
        !
        status = DftiComputeBackward(hand, ptr_snd)
        if(status /= 0) then
           write(*,*) "stopped in DftiComputeBackward", status
           stop
        endif
        ptr_snda( 1 : dffts%nr3x * dffts%nsw( dffts%mype + 1) ) &
             = ptr_snd( 1 : dffts%nr3x * dffts%nsw( dffts%mype + 1) )
        call tg_cft3s_scatter( ptr_snd, dffts, ptr_snda, 2, dffts%have_task_groups )
        !
        ! -- transport the send-pointer to the card
        !DIR$ offload_transfer target(mic:offload_device) in(ptr_snd:length(nnr) alloc_if(.false.) free_if(.false.))
     end if
     if(ibnd.gt.1) then
        !
        ! -- recover the receive buffer from the card
        !DIR$ offload_transfer target(mic:offload_device) out(ptr_rec:length(nnr) alloc_if(.false.) free_if(.false.))
        call tg_cft3s_scatter( ptr_rec, dffts,ptr_reca, -2, dffts%have_task_groups )
        tscale = 1.0_DP / ( dffts%nr3 )
        status = DftiComputeForward(hand, ptr_reca)
        if(status /= 0) then 
           write(*,*) "stopped in DftiComputeForward", status
           stop
        endif
        ptr_rec( 1 : dffts%nr3x * dffts%nsw( dffts%mype + 1) ) &
             =  ptr_reca( 1 : dffts%nr3x * dffts%nsw( dffts%mype + 1) )
        !$omp parallel do
        DO j = 1, n
           hpsi (j, ibnd-1)   = hpsi (j,ibnd-1)   + ptr_rec (nls(igk(j)))
        ENDDO
        !$omp end parallel do
     end if
     !
     ! -- if a signal has been issued (cahtwait=.true.) we need to wait for the signal
     t2=t2-dclock()
     if(catchwait) then
        !dir$ offload_wait target(mic:offload_device) wait(signal)
        catchwait=.false.
     end if
     t2=t2+dclock()
     !
     ! -- when all cycles have completed, we need to swap the buffers
     ptr_tmp=>ptr_rec
     ptr_rec=>ptr_cmp
     ptr_cmp=>ptr_snd
     ptr_snd=>ptr_tmp
     !
     ptr_tmp=>ptr_reca
     ptr_reca=>ptr_cmpa
     ptr_cmpa=>ptr_snda
     ptr_snda=>ptr_tmp
     !
  ENDDO

  ! -- deallocate all fields. start with the card allocations 
  !dir$ offload_transfer target(mic:offload_device) nocopy(psic_aux_1:length(nnr) alloc_if(.false.) free_if(.true.))
  !dir$ offload_transfer target(mic:offload_device) nocopy(psic_aux_2:length(nnr) alloc_if(.false.) free_if(.true.))
  !dir$ offload_transfer target(mic:offload_device) nocopy(psic_aux_3:length(nnr) alloc_if(.false.) free_if(.true.))
  deallocate(psic_aux_1)
  deallocate(aux_1)
  deallocate(psic_aux_2)
  deallocate(aux_2)
  deallocate(psic_aux_3)
  deallocate(aux_3)
  !
  ! -- free the pseudopotential field
  !dir$ offload_transfer target(mic:offload_device) nocopy(v:length(nnr) alloc_if(.false.) free_if(.true.))
  write(*,*) "FFT Time = ", dclock()-t1, t2
  RETURN
!END SUBROUTINE vloc_psi_k_catchwait
end SUBROUTINE vloc_psi_k
!
!


!-----------------------------------------------------------------------
SUBROUTINE vloc_psi_k_old(lda, n, m, psi, v, hpsi)
!SUBROUTINE vloc_psi_k(lda, n, m, psi, v, hpsi)
  !-----------------------------------------------------------------------
  !
  ! Calculation of Vloc*psi using dual-space technique - k-points
  !
  USE parallel_include
  USE kinds,   ONLY : DP
  USE gvecs, ONLY : nls, nlsm
  USE wvfct,   ONLY : igk
  USE mp_bands,      ONLY : me_bgrp
  USE fft_base,      ONLY : dffts, tg_gather
  USE fft_interfaces,ONLY : fwfft, invfft
  USE wavefunctions_module,  ONLY: psic
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: lda, n, m
  COMPLEX(DP), INTENT(in)   :: psi (lda, m)
  COMPLEX(DP), INTENT(inout):: hpsi (lda, m)
  REAL(DP), INTENT(in) :: v(dffts%nnr)
  !
  INTEGER :: ibnd, j, incr
  !
  LOGICAL :: use_tg
  ! Task Groups
  REAL(DP),    ALLOCATABLE :: tg_v(:)
  COMPLEX(DP), ALLOCATABLE :: tg_psic(:)
  INTEGER :: v_siz, idx, ioff
  !
#if defined(__CUDA) && !defined(__DISABLE_CUDA_VLOCPSI) && ( !defined(__PARA) || defined(__USE_3D_FFT) )
  CALL vloc_psi_k_gpu ( lda, n, m, psi, v, hpsi )
  RETURN
#endif
  !
  ! The following is dirty trick to prevent usage of task groups if
  ! the number of bands is smaller than the number of task groups 
  ! 
  use_tg = dffts%have_task_groups
  dffts%have_task_groups  = dffts%have_task_groups .and. ( m >= dffts%nogrp )
  !
  incr = 1
  !
  IF( dffts%have_task_groups ) THEN
     !
     v_siz =  dffts%tg_nnr * dffts%nogrp
     !
     ALLOCATE( tg_v   ( v_siz ) )
     ALLOCATE( tg_psic( v_siz ) )
     !
     CALL tg_gather( dffts, v, tg_v )
     incr = dffts%nogrp
     !
  ENDIF
  !
  ! the local potential V_Loc psi. First bring psi to real space
  !
  DO ibnd = 1, m, incr
     !
     IF( dffts%have_task_groups ) THEN
        !
        tg_psic = (0.d0, 0.d0)
        ioff   = 0
        !
        DO idx = 1, dffts%nogrp

           IF( idx + ibnd - 1 <= m ) THEN
!$omp parallel do
              DO j = 1, n
                 tg_psic(nls (igk(j))+ioff) =  psi(j,idx+ibnd-1)
              ENDDO
!$omp end parallel do
           ENDIF

           ioff = ioff + dffts%tg_nnr

        ENDDO
        !
        CALL  invfft ('Wave', tg_psic, dffts)
        !
     ELSE
        !
        psic(:) = (0.d0, 0.d0)
        psic (nls (igk(1:n))) = psi(1:n, ibnd)
        !
        CALL invfft ('Wave', psic, dffts)
        !
     ENDIF
     !
     !   fft to real space
     !   product with the potential v on the smooth grid
     !   back to reciprocal space
     !
     IF( dffts%have_task_groups ) THEN
        !
!$omp parallel do
        DO j = 1, dffts%nr1x*dffts%nr2x*dffts%tg_npp( me_bgrp + 1 )
           tg_psic (j) = tg_psic (j) * tg_v(j)
        ENDDO
!$omp end parallel do
        !
        CALL fwfft ('Wave',  tg_psic, dffts)
        !
     ELSE
        !
!$omp parallel do
        DO j = 1, dffts%nnr
           psic (j) = psic (j) * v(j)
        ENDDO
!$omp end parallel do
        !
        CALL fwfft ('Wave', psic, dffts)
        !
     ENDIF
     !
     !   addition to the total product
     !
     IF( dffts%have_task_groups ) THEN
        !
        ioff   = 0
        !
        DO idx = 1, dffts%nogrp
           !
           IF( idx + ibnd - 1 <= m ) THEN
!$omp parallel do
              DO j = 1, n
                 hpsi (j, ibnd+idx-1) = hpsi (j, ibnd+idx-1) + tg_psic( nls(igk(j)) + ioff )
              ENDDO
!$omp end parallel do
           ENDIF
           !
           ioff = ioff + dffts%nr3x * dffts%nsw( me_bgrp + 1 )
           !
        ENDDO
        !
     ELSE
!$omp parallel do
        DO j = 1, n
           hpsi (j, ibnd)   = hpsi (j, ibnd)   + psic (nls(igk(j)))
        ENDDO
!$omp end parallel do
     ENDIF
     !
  ENDDO
  !
  IF( dffts%have_task_groups ) THEN
     !
     DEALLOCATE( tg_psic )
     DEALLOCATE( tg_v )
     !
  ENDIF
  dffts%have_task_groups = use_tg
  !
  RETURN
END SUBROUTINE VLOC_PSI_K_OLD
!end SUBROUTINE vloc_psi_k
!
!-----------------------------------------------------------------------
SUBROUTINE vloc_psi_nc (lda, n, m, psi, v, hpsi)
  !-----------------------------------------------------------------------
  !
  ! Calculation of Vloc*psi using dual-space technique - noncolinear
  !
  USE parallel_include
  USE kinds,   ONLY : DP
  USE gvecs, ONLY : nls, nlsm
  USE wvfct,   ONLY : igk
  USE mp_bands,      ONLY : me_bgrp
  USE fft_base,      ONLY : dffts, dfftp, tg_gather
  USE fft_interfaces,ONLY : fwfft, invfft
  USE lsda_mod,      ONLY : nspin
  USE spin_orb,      ONLY : domag
  USE noncollin_module,     ONLY: npol
  USE wavefunctions_module, ONLY: psic_nc
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: lda, n, m
  REAL(DP), INTENT(in) :: v(dfftp%nnr,4) ! beware dimensions!
  COMPLEX(DP), INTENT(in)   :: psi (lda*npol, m)
  COMPLEX(DP), INTENT(inout):: hpsi (lda,npol,m)
  !
  INTEGER :: ibnd, j,ipol, incr, is
  COMPLEX(DP) :: sup, sdwn
  !
  LOGICAL :: use_tg
  ! Variables for task groups
  REAL(DP),    ALLOCATABLE :: tg_v(:,:)
  COMPLEX(DP), ALLOCATABLE :: tg_psic(:,:)
  INTEGER :: v_siz, idx, ioff
  !
  !
  incr = 1
  !
  ! The following is dirty trick to prevent usage of task groups if
  ! the number of bands is smaller than the number of task groups 
  ! 
  use_tg = dffts%have_task_groups
  dffts%have_task_groups  = dffts%have_task_groups .and. ( m >= dffts%nogrp )
  !
  IF( dffts%have_task_groups ) THEN
     v_siz = dffts%tg_nnr * dffts%nogrp
     IF (domag) THEN
        ALLOCATE( tg_v( v_siz, 4 ) )
        DO is=1,nspin
           CALL tg_gather( dffts, v(:,is), tg_v(:,is) )
        ENDDO
     ELSE
        ALLOCATE( tg_v( v_siz, 1 ) )
        CALL tg_gather( dffts, v(:,1), tg_v(:,1) )
     ENDIF
     ALLOCATE( tg_psic( v_siz, npol ) )
     incr = dffts%nogrp
  ENDIF
  !
  ! the local potential V_Loc psi. First the psi in real space
  !
  DO ibnd = 1, m, incr

     IF( dffts%have_task_groups ) THEN
        !
        DO ipol = 1, npol
           !
           tg_psic(:,ipol) = ( 0.D0, 0.D0 )
           ioff   = 0
           !
           DO idx = 1, dffts%nogrp
              !
              IF( idx + ibnd - 1 <= m ) THEN
                 DO j = 1, n
                    tg_psic( nls( igk(j) ) + ioff, ipol ) = psi( j +(ipol-1)*lda, idx+ibnd-1 )
                 ENDDO
              ENDIF

              ioff = ioff + dffts%tg_nnr

           ENDDO
           !
           CALL invfft ('Wave', tg_psic(:,ipol), dffts)
           !
        ENDDO
        !
     ELSE
        psic_nc = (0.d0,0.d0)
        DO ipol=1,npol
           DO j = 1, n
              psic_nc(nls(igk(j)),ipol) = psi(j+(ipol-1)*lda,ibnd)
           ENDDO
           CALL invfft ('Wave', psic_nc(:,ipol), dffts)
        ENDDO
     ENDIF

     !
     !   product with the potential v = (vltot+vr) on the smooth grid
     !
     IF( dffts%have_task_groups ) THEN
        IF (domag) THEN
           DO j=1, dffts%nr1x*dffts%nr2x*dffts%tg_npp( me_bgrp + 1 )
              sup = tg_psic(j,1) * (tg_v(j,1)+tg_v(j,4)) + &
                    tg_psic(j,2) * (tg_v(j,2)-(0.d0,1.d0)*tg_v(j,3))
              sdwn = tg_psic(j,2) * (tg_v(j,1)-tg_v(j,4)) + &
                     tg_psic(j,1) * (tg_v(j,2)+(0.d0,1.d0)*tg_v(j,3))
              tg_psic(j,1)=sup
              tg_psic(j,2)=sdwn
           ENDDO
        ELSE
           DO j=1, dffts%nr1x*dffts%nr2x*dffts%tg_npp( me_bgrp + 1 )
              tg_psic(j,:) = tg_psic(j,:) * tg_v(j,1)
           ENDDO
        ENDIF
     ELSE
        IF (domag) THEN
           DO j=1, dffts%nnr
              sup = psic_nc(j,1) * (v(j,1)+v(j,4)) + &
                    psic_nc(j,2) * (v(j,2)-(0.d0,1.d0)*v(j,3))
              sdwn = psic_nc(j,2) * (v(j,1)-v(j,4)) + &
                     psic_nc(j,1) * (v(j,2)+(0.d0,1.d0)*v(j,3))
              psic_nc(j,1)=sup
              psic_nc(j,2)=sdwn
           ENDDO
        ELSE
           DO j=1, dffts%nnr
              psic_nc(j,:) = psic_nc(j,:) * v(j,1)
           ENDDO
        ENDIF
     ENDIF
     !
     !   back to reciprocal space
     !
     IF( dffts%have_task_groups ) THEN
        !
        DO ipol = 1, npol

           CALL fwfft ('Wave', tg_psic(:,ipol), dffts)
           !
           ioff   = 0
           !
           DO idx = 1, dffts%nogrp
              !
              IF( idx + ibnd - 1 <= m ) THEN
                 DO j = 1, n
                    hpsi (j, ipol, ibnd+idx-1) = hpsi (j, ipol, ibnd+idx-1) + &
                                           tg_psic( nls(igk(j)) + ioff, ipol )
                 ENDDO
              ENDIF
              !
              ioff = ioff + dffts%nr3x * dffts%nsw( me_bgrp + 1 )
              !
           ENDDO

        ENDDO
        !
     ELSE

        DO ipol=1,npol
           CALL fwfft ('Wave', psic_nc(:,ipol), dffts)
        ENDDO
        !
        !   addition to the total product
        !
        DO ipol=1,npol
           DO j = 1, n
              hpsi(j,ipol,ibnd) = hpsi(j,ipol,ibnd) + psic_nc(nls(igk(j)),ipol)
           ENDDO
        ENDDO

     ENDIF

  ENDDO

  IF( dffts%have_task_groups ) THEN
     !
     DEALLOCATE( tg_v )
     DEALLOCATE( tg_psic )
     !
  ENDIF
  dffts%have_task_groups = use_tg
  !
  RETURN
END SUBROUTINE vloc_psi_nc
