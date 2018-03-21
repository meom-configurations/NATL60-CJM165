MODULE trcdiasub
   !!======================================================================
   !!                     ***  MODULE  trcdiasub  ***
   !! TOP :  computes passive tracer subduction 
   !!=====================================================================
   !! History :   1.0  !  2006-06  (P. Karleskind)  original code
   !!              -   !  2011-10  (P. Karleskind)  F90 
   !!----------------------------------------------------------------------
#if  defined key_top  &&  defined key_diasub
   !!----------------------------------------------------------------------
   !!   'key_top'  and  'key_diasub'        TOP model + passive tracer sub
   !!----------------------------------------------------------------------
   !!   trc_sub      : computes passive tracer subdcution
   !!----------------------------------------------------------------------
   USE oce_trc
   USE par_trc
   USE trc
   USE lib_print
   USE iom
   USE diasub

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_dia_sub_adv
   PUBLIC   trc_dia_sub_adv_eiv
   PUBLIC   trc_dia_sub_ldf
   PUBLIC   trc_dia_sub_zdf
   PUBLIC   trc_dia_sub_alloc

   REAL(wp), PUBLIC, SAVE ,DIMENSION(:,:,:,:), ALLOCATABLE ::   trsub  !: tracer subduction array

   !! * Substitutions
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcsub.F90 2715 2011-03-30 15:58:35Z rblod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION trc_dia_sub_alloc()
      !!---------------------------------------------------------------------
      !
      ALLOCATE( trsub(jpi,jpj,jptra,jptrsub ), STAT=trc_dia_sub_alloc)
      !
      IF( lk_mpp           )   CALL mpp_sum ( trc_dia_sub_alloc )
      IF( trc_dia_sub_alloc /= 0)   CALL ctl_warn('trc_dia_sub_alloc: failed to allocate arrays.')
      !
   END FUNCTION trc_dia_sub_alloc

   SUBROUTINE trc_dia_sub_adv(kt, pun, pvn, pwn ) 
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE trc_sub  ***
      !!
      !! ** Purpose :   computes passive tracer subdcution for mld and adv (without eiv) 
      !!
      !! ** Method  :   
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt     ! ocean time-step
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in) ::   pun, pvn, pwn   ! 3 ocean velocity components 
      !!
      INTEGER  ::  ji, jj, jn, jk 
      INTEGER  ::  ikn, ikni, iknj, ikw, iku, ikv, ikb  
      REAL(wp) ::  zuints, zvints, zhint, ztrsu, ztrsv
      !!----------------------------------------------------------------------

!     trsub(:,:,:,:)=0.


      DO jn = 1,jptra

! Add tracers concentration
         DO jj = fs_2, jpj
            DO ji = fs_2, jpi
               ikn   = nmln(ji,jj)  ;  ikni = nmln(ji-1,jj)  ;  iknj = nmln(ji,jj-1)
               ikw   = nkind(ji,jj) 
! W term
               trsub(ji,jj,jn,jpsub_zad) = -pwn(ji,jj,ikn) * trn(ji,jj,ikw,jn) * tmask(ji,jj,ikw)
! Adv term
               DO jk = minku(ji,jj), maxku(ji,jj)
                   iku  = nindui(ji,jj,jk) 
                   ztrsu = umask(ji-1,jj,jk)      &
                   &    * sign( 1, ikni - ikn )  &
                   &    * pun(ji-1,jj,jk) * trn(iku,jj,jk,jn) * tmask(iku,jj,jk)                     
                  trsub(ji,jj,jn,jpsub_xad) =  trsub(ji,jj,jn,jpsub_xad) + ztrsu 
               ENDDO
               !
               DO jk = minkv(ji,jj), maxkv(ji,jj)
                  ikv = nindvj(ji,jj,jk) 
                  ztrsv = vmask(ji,jj-1,jk)      &
                   &    * sign( 1, iknj - ikn )  &
                   &    * pvn(ji,jj-1,jk) * trn(ji,ikv,jk,jn) * tmask(ji,ikv,jk)                     
                  trsub(ji,jj,jn,jpsub_yad) =  trsub(ji,jj,jn,jpsub_yad) + ztrsv 
               ENDDO
               !
              IF( kt > nit000 ) THEN
                 ikn = nmln(ji,jj)   ;     ikb = nmlb(ji,jj) 
                 zhint = 0._wp
                 IF( ikn /= ikb ) THEN
                   DO jk = minh(ji,jj), maxh(ji,jj)
                      zhint = zhint +  gmlh(ji,jj,jk) * trb(ji,jj,jk,jn)
                   END DO
                   trsub(ji,jj,jn,jpsub_mld) = -zhint
                 ENDIF
             ENDIF

            END DO
         END DO
      END DO

      hmlpb(:,:) = hmlp(:,:)
      nmlb (:,:) = nmln(:,:)
      !
   END SUBROUTINE trc_dia_sub_adv


   SUBROUTINE trc_dia_sub_adv_eiv( kt, pun, pvn, pwn )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE trc_dia_sub_adv_eiv  ***
      !!
      !! ** Purpose :   computes passive tracer subdcution for adv+eiv
      !!
      !! ** Method  :
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt     ! ocean time-step
      REAL(wp), DIMENSION(jpi,jpj,jpk     ), INTENT(in   ) ::   pun, pvn, pwn   ! 3 ocean velocity components
#if defined key_trcldf_eiv
      !!
      INTEGER ::   ji, jj, jn, jk
      INTEGER ::   ikn, ikb, ikni, iknj
      INTEGER ::   iieiu,iieiv, iieiw
      REAL(wp) ::  ztrsu, ztrsv
      INTEGER, POINTER, DIMENSION(:,:)  :: ikindeiu, ikindeiv, ikindeiw
      !!----------------------------------------------------------------------
      !
      IF( kt == nit000 )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'trc_dia_sub_adv_eiv : passive tracer subduction for eiv '
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~ '
      ENDIF
      !

  !   trsub(:,:,:,:)=0.


      CALL wrk_alloc(jpi,jpj, ikindeiu, ikindeiv, ikindeiw )

      DO jj = fs_2, jpj
         DO ji = fs_2, jpi
            ikn  = nmln(ji,jj)  ;  ikni = nmln(ji-1,jj)  ;  iknj = nmln(ji,jj-1)
            !
            IF( pwn(ji,jj,ikn) < 0._wp ) THEN     ;  ikindeiw(ji,jj) = ikn - 1
            ELSE                                  ;  ikindeiw(ji,jj) = ikn
            ENDIF
            !
            DO jk = minku(ji,jj), maxku(ji,jj)
              IF( pun(ji-1,jj,jk) < 0._wp ) THEN  ;  ikindeiu(ji,jj) = ji
              ELSE                                ;  ikindeiu(ji,jj) = ji - 1
              ENDIF
            ENDDO
            !
            DO jk = minkv(ji,jj), maxkv(ji,jj)
              IF( pvn(ji,jj-1,jk) < 0._wp ) THEN  ;  ikindeiv(ji,jj) = jj
              ELSE                                ;  ikindeiv(ji,jj) = jj - 1
             ENDIF
            ENDDO
     !
            END DO
         END DO


      DO jn = 1,jptra

! Add tracers concentration
         DO jj = fs_2, jpj
            DO ji = fs_2, jpi
               ikn   = nmln(ji,jj)  ;  ikni = nmln(ji-1,jj)  ;  iknj = nmln(ji,jj-1)
               iieiu = ikindeiu(ji,jj) 
               iieiv = ikindeiv(ji,jj) 
               iieiw = ikindeiw(ji,jj) 
! W term
               trsub(ji,jj,jn,jpsub_zei) = -pwn(ji,jj,ikn) * trn(ji,jj,iieiw,jn) * tmask(ji,jj,ikn)
! Adv term
               DO jk = minku(ji,jj), maxku(ji,jj)
                  ztrsu = umask(ji-1,jj,jk)      &
                   &    * sign( 1, ikni - ikn )  &
                   &    * pun(ji-1,jj,jk) * trn(iieiu,jj,jk,jn) * tmask(iieiu,jj,jk)                     
                  trsub(ji,jj,jn,jpsub_xei) =  trsub(ji,jj,jn,jpsub_xei) + ztrsu 
               ENDDO
               !
               DO jk = minkv(ji,jj), maxkv(ji,jj)
                  ztrsv = vmask(ji,jj-1,jk)      &
                   &    * sign( 1, iknj - ikn )  &
                   &    * pvn(ji,jj-1,jk) * trn(ji,iieiv,jk,jn) * tmask(ji,iieiv,jk)                     
                  trsub(ji,jj,jn,jpsub_yei) =  trsub(ji,jj,jn,jpsub_yei) + ztrsv 
               ENDDO
               !
            END DO
         END DO
      END DO
      CALL wrk_dealloc(jpi,jpj, ikindeiu, ikindeiv, ikindeiw )
#endif
      !
   END SUBROUTINE trc_dia_sub_adv_eiv


   SUBROUTINE trc_dia_sub_ldf( kt, ptru, ptrv, ptrw )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE trc_dia_sub_ldf  ***
      !!
      !! ** Purpose :   computes passive tracer subduction for ldf
      !!
      !! ** Method  :
      !!----------------------------------------------------------------------
      INTEGER                              , INTENT(in   ) ::   kt         ! ocean time-step index
      REAL(wp), DIMENSION(jpi,jpj,jpk,jptra), INTENT(in   ) ::   ptru, ptrv, ptrw   ! 3 ocean lateral diffusion trends components
      !!
      INTEGER            :: ji, jj, jk, jn
      INTEGER            :: ikn, ikni, iknj, iknp1
      !!----------------------------------------------------------------------
      !
      IF( kt == nit000 )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'trc_dia_sub_ldf : passive tracer subduction for lateral diffusion '
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~ '
      ENDIF
      !
      DO jn = 1, jptra
         DO ji = fs_2, jpi
            DO jj = fs_2, jpj
               ikn  = nmln(ji,jj)  ;  ikni = nmln(ji-1,jj)  ;  iknj = nmln(ji,jj-1)
               !
               DO jk = minku(ji,jj),  maxku(ji,jj)
                  trsub(ji,jj,jn,jpsub_xlf)  = trsub(ji,jj,jn,jpsub_xlf) &
                   &                         + umask(ji-1,jj,jk) * SIGN( 1, ikni - ikn ) * ptru(ji-1,jj,jk,jn)
               END DO
               !
               DO jk = minkv(ji,jj), maxkv(ji,jj)
                   trsub(ji,jj,jn,jpsub_ylf) =  trsub(ji,jj,jn,jpsub_ylf) &
                    &                        + vmask(ji,jj-1,jk) * SIGN( 1, iknj - ikn ) * ptrv(ji,jj-1,jk,jn)
               END DO
            END DO
         END DO
         !
         DO jj = 1, jpj
            DO ji = 1, jpi
               iknp1 = nmln(ji,jj) + 1
               trsub(ji,jj,jn,jpsub_zlf) = ptrw(ji,jj,iknp1,jn)
            END DO
         END DO
      END DO
      !
   END SUBROUTINE trc_dia_sub_ldf

   SUBROUTINE trc_dia_sub_zdf( kt, ptrzdf )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE trc_dia_sub_zdf  ***
      !!
      !! ** Purpose :   computes passive tracer subduction for zdf
      !!
      !! ** Method  :
      !!----------------------------------------------------------------------
      INTEGER                              , INTENT(in   ) ::   kt         ! ocean time-step index
      REAL(wp), DIMENSION(jpi,jpj,jpk,jptra), INTENT(in) :: ptrzdf   !  ocean vertical diffusion trend
      INTEGER            :: ji, jj, jk, jn
      !!----------------------------------------------------------------------
      !
      IF( kt == nit000 )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'trc_dia_sub_zdf : passive tracer subduction for vertical diffusion '
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~ '
      ENDIF
      !
      DO jn = 1, jptra
         DO ji = 2, jpi
            DO jj = 2, jpj
               DO jk = 1, nkind(ji,jj)
                  trsub(ji,jj,jn,jpsub_zdf) = trsub(ji,jj,jn,jpsub_zdf)  &
                      &                     - fse3t(ji,jj,jk) * e2t(ji,jj) * e1t(ji,jj) * ptrzdf(ji,jj,jk,jn)

               END DO
            END DO
         END DO
      END DO
      !
   END SUBROUTINE trc_dia_sub_zdf

#else
   !!----------------------------------------------------------------------
   !!   Dummy module                              NO passive tracer Subduction
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_dia_sub_adv( kt, pun, pvn, pwn )              ! Empty routine
      INTEGER  ::   kt
      REAL, DIMENSION(:,:,:) ::   pun, pvn, pwn
      WRITE(*,*) 'trc_dia_sub_adv: You should not have seen this print! error?', kt, pun(1,1,1), pvn(1,1,1), pwn(1,1,1)
   END SUBROUTINE trc_dia_sub_adv

   SUBROUTINE trc_dia_sub_adv_eiv( kt, pun, pvn, pwn )              ! Empty routine
      INTEGER  ::   kt
      REAL, DIMENSION(:,:,:) ::   pun, pvn, pwn
      WRITE(*,*) 'trc_dia_sub_adv_eiv: You should not have seen this print! error?', kt, pun(1,1,1), pvn(1,1,1), pwn(1,1,1)
   END SUBROUTINE trc_dia_sub_adv_eiv

   SUBROUTINE trc_dia_sub_ldf( kt,  ptru, ptrv, ptrw )              ! Empty routine
      REAL, DIMENSION(:,:,:,:) ::   ptru, ptrv, ptrw
      WRITE(*,*) 'trc_dia_sub_ldf: You should not have seen this print! error?', ptru(1,1,1,1), ptrv(1,1,1,1), ptrw(1,1,1,1)
   END SUBROUTINE trc_dia_sub_ldf

   SUBROUTINE trc_dia_sub_zdf( kt,  ptrzdf )              ! Empty routine
      REAL, DIMENSION(:,:,:,:) ::   ptrzdf
      WRITE(*,*) 'trc_dia_sub_zdf: You should not have seen this print! error?', ptrzdf(1,1,1,1)
   END SUBROUTINE trc_dia_sub_zdf

#endif

   !!======================================================================
END MODULE trcdiasub
