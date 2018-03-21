MODULE traldf_lap
   !!==============================================================================
   !!                       ***  MODULE  traldf_lap  ***
   !! Ocean  tracers:  horizontal component of the lateral tracer mixing trend
   !!==============================================================================
   !! History :  OPA  !  87-06  (P. Andrich, D. L Hostis)  Original code
   !!                 !  91-11  (G. Madec)
   !!                 !  95-11  (G. Madec)  suppress volumetric scale factors
   !!                 !  96-01  (G. Madec)  statement function for e3
   !!            NEMO !  02-06  (G. Madec)  F90: Free form and module
   !!            1.0  !  04-08  (C. Talandier) New trends organization
   !!                 !  05-11  (G. Madec)  add zps case
   !!            3.0  !  10-06  (C. Ethe, G. Madec) Merge TRA-TRC
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   tra_ldf_lap  : update the tracer trend with the horizontal diffusion
   !!                 using a iso-level harmonic (laplacien) operator.
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and active tracers
   USE dom_oce         ! ocean space and time domain
   USE ldftra_oce      ! ocean active tracers: lateral physics
   USE in_out_manager  ! I/O manager
   USE diaptr          ! poleward transport diagnostics
   USE trc_oce         ! share passive tracers/Ocean variables
   USE lib_mpp         ! MPP library
   USE timing          ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tra_ldf_lap   ! routine called by step.F90

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "ldftra_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: traldf_lap.F90 5147 2015-03-13 10:01:32Z cetlod $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE tra_ldf_lap( kt, kit000, cdtype, pgu , pgv ,    &
      &                                        pgui, pgvi,    &
      &                                ptb, pta, kjpt,        &
      &                                ptru, ptrv, ptrw ) 
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_ldf_lap  ***
      !!                   
      !! ** Purpose :   Compute the before horizontal tracer (t & s) diffusive 
      !!      trend and add it to the general trend of tracer equation.
      !!
      !! ** Method  :   Second order diffusive operator evaluated using before
      !!      fields (forward time scheme). The horizontal diffusive trends of 
      !!      the tracer is given by:
      !!          difft = 1/(e1t*e2t*e3t) {  di-1[ aht e2u*e3u/e1u di(tb) ]
      !!                                   + dj-1[ aht e1v*e3v/e2v dj(tb) ] }
      !!      Add this trend to the general tracer trend pta :
      !!          pta = pta + difft
      !!
      !! ** Action  : - Update pta arrays with the before iso-level 
      !!                harmonic mixing trend.
      !!----------------------------------------------------------------------
      USE oce, ONLY:   ztu => ua , ztv => va  ! (ua,va) used as workspace
      !
      INTEGER                              , INTENT(in   ) ::   kt         ! ocean time-step index
      INTEGER                              , INTENT(in   ) ::   kit000          ! first time step index
      CHARACTER(len=3)                     , INTENT(in   ) ::   cdtype     ! =TRA or TRC (tracer indicator)
      INTEGER                              , INTENT(in   ) ::   kjpt       ! number of tracers
      REAL(wp), DIMENSION(jpi,jpj    ,kjpt), INTENT(in   ) ::   pgu, pgv   ! tracer gradient at pstep levels
      REAL(wp), DIMENSION(jpi,jpj,    kjpt), INTENT(in   ) ::   pgui, pgvi ! tracer gradient at top levels
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(in   ) ::   ptb        ! before and now tracer fields
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(inout) ::   pta        ! tracer trend 
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(out  ), OPTIONAL ::   ptru, ptrv, ptrw ! optional u/v/w-tracer trend 

      !
      INTEGER  ::   ji, jj, jk, jn       ! dummy loop indices
      INTEGER  ::   iku, ikv, ierr       ! local integers
      REAL(wp) ::   zabe1, zabe2, zbtr   ! local scalars
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 ) CALL timing_start('tra_ldf_lap')
      !
      IF( kt == kit000 )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_ldf_lap : iso-level laplacian diffusion on ', cdtype
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~ '
      ENDIF

      !                                                          ! =========== !
      DO jn = 1, kjpt                                            ! tracer loop !
         !                                                       ! =========== !    
         DO jk = 1, jpkm1                                            ! slab loop
            !                                           
            ! 1. First derivative (gradient)
            ! -------------------
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1   ! vector opt.
                  zabe1 = fsahtu(ji,jj,jk) * umask(ji,jj,jk) * re2u_e1u(ji,jj) * fse3u_n(ji,jj,jk)
                  zabe2 = fsahtv(ji,jj,jk) * vmask(ji,jj,jk) * re1v_e2v(ji,jj) * fse3v_n(ji,jj,jk)
                  ztu(ji,jj,jk) = zabe1 * ( ptb(ji+1,jj  ,jk,jn) - ptb(ji,jj,jk,jn) )
                  ztv(ji,jj,jk) = zabe2 * ( ptb(ji  ,jj+1,jk,jn) - ptb(ji,jj,jk,jn) )
               END DO
            END DO
            IF( ln_zps ) THEN      ! set gradient at partial step level for the last ocean cell
               DO jj = 1, jpjm1
                  DO ji = 1, fs_jpim1   ! vector opt.
                     ! last level
                     iku = mbku(ji,jj)
                     ikv = mbkv(ji,jj)
                     IF( iku == jk ) THEN
                        zabe1 = fsahtu(ji,jj,iku) * umask(ji,jj,iku) * re2u_e1u(ji,jj) * fse3u_n(ji,jj,iku)
                        ztu(ji,jj,jk) = zabe1 * pgu(ji,jj,jn)
                     ENDIF
                     IF( ikv == jk ) THEN
                        zabe2 = fsahtv(ji,jj,ikv) * vmask(ji,jj,ikv) * re1v_e2v(ji,jj) * fse3v_n(ji,jj,ikv)
                        ztv(ji,jj,jk) = zabe2 * pgv(ji,jj,jn)
                     ENDIF
                  END DO
               END DO
            ENDIF
            ! (ISH)
            IF( ln_zps .AND. ln_isfcav ) THEN      ! set gradient at partial step level for the first ocean cell
                                                   ! into a cavity
               DO jj = 1, jpjm1
                  DO ji = 1, fs_jpim1   ! vector opt.
                     ! ice shelf level level MAX(2,jk) => only where ice shelf
                     iku = miku(ji,jj) 
                     ikv = mikv(ji,jj) 
                     IF( iku == MAX(2,jk) ) THEN 
                        zabe1 = fsahtu(ji,jj,iku) * umask(ji,jj,iku) * re2u_e1u(ji,jj) * fse3u_n(ji,jj,iku) 
                        ztu(ji,jj,jk) = zabe1 * pgui(ji,jj,jn) 
                     ENDIF 
                     IF( ikv == MAX(2,jk) ) THEN 
                        zabe2 = fsahtv(ji,jj,ikv) * vmask(ji,jj,ikv) * re1v_e2v(ji,jj) * fse3v_n(ji,jj,ikv) 
                        ztv(ji,jj,jk) = zabe2 * pgvi(ji,jj,jn) 
                     END IF 
                  END DO
               END DO
            ENDIF
         
         
            ! 2. Second derivative (divergence) added to the general tracer trends
            ! ---------------------------------------------------------------------
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zbtr = 1._wp / ( e12t(ji,jj) * fse3t_n(ji,jj,jk) )
                  ! horizontal diffusive trends added to the general tracer trends
                  pta(ji,jj,jk,jn) = pta(ji,jj,jk,jn) + zbtr * (  ztu(ji,jj,jk) - ztu(ji-1,jj,jk)   &
                     &                                          + ztv(ji,jj,jk) - ztv(ji,jj-1,jk)  )
               END DO
            END DO
            !
         END DO                                             !  End of slab  
         !
         ! "Poleward" diffusive heat or salt transports
         IF( cdtype == 'TRA' .AND. ln_diaptr ) THEN
            IF( jn  == jp_tem)   htr_ldf(:) = ptr_sj( ztv(:,:,:) )
            IF( jn  == jp_sal)   str_ldf(:) = ptr_sj( ztv(:,:,:) )
         ENDIF
         !   
         IF( PRESENT( ptru ) ) THEN
            ptru(:,:,:,jn) = ztu(:,:,:)
            ptrv(:,:,:,jn) = ztv(:,:,:)
            ptrw(:,:,:,jn) = 0.
         ENDIF
         !
         !                                                  ! ==================
      END DO                                                ! end of tracer loop
      !                                                     ! ==================
      IF( nn_timing == 1 ) CALL timing_stop('tra_ldf_lap')
      !
   END SUBROUTINE tra_ldf_lap

   !!==============================================================================
END MODULE traldf_lap
