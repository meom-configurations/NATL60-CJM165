MODULE traldf_bilap
   !!==============================================================================
   !!                   ***  MODULE  traldf_bilap  ***
   !! Ocean  tracers:  horizontal component of the lateral tracer mixing trend
   !!==============================================================================
   !! History :  OPA  !  1991-11  (G. Madec)  Original code
   !!                 !  1993-03  (M. Guyon)  symetrical conditions
   !!                 !  1995-11  (G. Madec)  suppress volumetric scale factors
   !!                 !  1996-01  (G. Madec)  statement function for e3
   !!                 !  1996-01  (M. Imbard)  mpp exchange
   !!                 !  1997-07  (G. Madec)  optimization, and ahtt
   !!            8.5  !  2002-08  (G. Madec)  F90: Free form and module
   !!   NEMO     1.0  !  2004-08  (C. Talandier) New trends organization
   !!             -   !  2005-11  (G. Madec)  zps or sco as default option
   !!            3.3  !  2010-05  (C. Ethe, G. Madec)  merge TRC-TRA 
   !!==============================================================================

   !!----------------------------------------------------------------------
   !!   tra_ldf_bilap : update the tracer trend with the horizontal diffusion
   !!                   using a iso-level biharmonic operator
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and active tracers
   USE dom_oce         ! ocean space and time domain
   USE ldftra_oce      ! ocean tracer   lateral physics
   USE in_out_manager  ! I/O manager
   USE ldfslp          ! iso-neutral slopes 
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE diaptr          ! poleward transport diagnostics
   USE trc_oce         ! share passive tracers/Ocean variables
   USE lib_mpp         ! MPP library
   USE wrk_nemo       ! Memory Allocation
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tra_ldf_bilap   ! routine called by step.F90

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "ldftra_substitute.h90"
#  include "ldfeiv_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: traldf_bilap.F90 5147 2015-03-13 10:01:32Z cetlod $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS
 
   SUBROUTINE tra_ldf_bilap( kt, kit000, cdtype, pgu, pgv,            &
      &                                          pgui, pgvi,          &
      &                                  ptb, pta, kjpt,              &
      &                                  ptru, ptrv, ptrw           )  
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_ldf_bilap  ***
      !!
      !! ** Purpose :   Compute the before horizontal tracer diffusive 
      !!      trend and add it to the general trend of tracer equation.
      !!
      !! ** Method  :   4th order diffusive operator along model level surfaces 
      !!      evaluated using before fields (forward time scheme). The hor.
      !!      diffusive trends  is given by:
      !!      Laplacian of tb:
      !!         zlt   = 1/(e1t*e2t*e3t) {  di-1[ e2u*e3u/e1u di(tb) ]
      !!                                  + dj-1[ e1v*e3v/e2v dj(tb) ]  }
      !!      Multiply by the eddy diffusivity coef. and insure lateral bc:
      !!        zlt   = ahtt * zlt
      !!        call to lbc_lnk
      !!      Bilaplacian (laplacian of zlt):
      !!         difft = 1/(e1t*e2t*e3t) {  di-1[ e2u*e3u/e1u di(zlt) ]
      !!                                  + dj-1[ e1v*e3v/e2v dj(zlt) ]  }
      !!
      !!      Add this trend to the general trend
      !!         (pta) = (pta) + ( difft )
      !!
      !! ** Action : - Update pta arrays with the before iso-level
      !!               biharmonic mixing trend.
      !!----------------------------------------------------------------------
      USE oce     , ONLY:   ztu  => ua       , ztv  => va                           ! (ua,va) used as workspace
      !!
      INTEGER                              , INTENT(in   ) ::   kt         ! ocean time-step index
      INTEGER                              , INTENT(in   ) ::   kit000          ! first time step index
      CHARACTER(len=3)                     , INTENT(in   ) ::   cdtype     ! =TRA or TRC (tracer indicator)
      INTEGER                              , INTENT(in   ) ::   kjpt       ! number of tracers
      REAL(wp), DIMENSION(jpi,jpj,    kjpt), INTENT(in   ) ::   pgu , pgv  ! tracer gradient at pstep levels
      REAL(wp), DIMENSION(jpi,jpj,    kjpt), INTENT(in   ) ::   pgui, pgvi ! tracer gradient at pstep levels
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(in   ) ::   ptb        ! before and now tracer fields
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(inout) ::   pta        ! tracer trend 
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(out  ), OPTIONAL ::   ptru, ptrv, ptrw ! optional u/v/w-tracer trend 

      !!
      INTEGER  ::  ji, jj, jk, jn   ! dummy loop indices
      REAL(wp) ::  zbtr, ztra       ! local scalars
      REAL(wp), POINTER, DIMENSION(:,:) ::  zeeu, zeev, zlt
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start( 'tra_ldf_bilap')
      !
      CALL wrk_alloc( jpi, jpj, zeeu, zeev, zlt ) 
      !

      IF( kt == kit000 )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_ldf_bilap : iso-level biharmonic operator on ', cdtype
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~'
      ENDIF
      !                                                          ! ===========
      DO jn = 1, kjpt                                            ! tracer loop
         !                                                       ! ===========
         !                                               
         DO jk = 1, jpkm1                                        ! Horizontal slab
            !                                             
            !                          !==  Initialization of metric arrays (for z- or s-coordinates)  ==!
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1   ! vector opt.
                  zeeu(ji,jj) = re2u_e1u(ji,jj) * fse3u_n(ji,jj,jk) * umask(ji,jj,jk)
                  zeev(ji,jj) = re1v_e2v(ji,jj) * fse3v_n(ji,jj,jk) * vmask(ji,jj,jk)
               END DO
            END DO
            !                          !==  Laplacian  ==!
            !
            DO jj = 1, jpjm1                 ! First derivative (gradient)
               DO ji = 1, fs_jpim1   ! vector opt.
                  ztu(ji,jj,jk) = zeeu(ji,jj) * ( ptb(ji+1,jj  ,jk,jn) - ptb(ji,jj,jk,jn) )
                  ztv(ji,jj,jk) = zeev(ji,jj) * ( ptb(ji  ,jj+1,jk,jn) - ptb(ji,jj,jk,jn) )
               END DO
            END DO
            !
            IF( ln_zps ) THEN                ! set gradient at partial step level (last ocean level)
               DO jj = 1, jpjm1
                  DO ji = 1, jpim1
                     IF( mbku(ji,jj) == jk )  ztu(ji,jj,jk) = zeeu(ji,jj) * pgu(ji,jj,jn)
                     IF( mbkv(ji,jj) == jk )  ztv(ji,jj,jk) = zeev(ji,jj) * pgv(ji,jj,jn)
                  END DO
               END DO
            ENDIF
            ! (ISH)
            IF( ln_zps .AND. ln_isfcav ) THEN ! set gradient at partial step level (first ocean level in a cavity)
               DO jj = 1, jpjm1
                  DO ji = 1, jpim1
                     IF( miku(ji,jj) == MAX(jk,2) )  ztu(ji,jj,jk) = zeeu(ji,jj) * pgui(ji,jj,jn)
                     IF( mikv(ji,jj) == MAX(jk,2) )  ztu(ji,jj,jk) = zeev(ji,jj) * pgvi(ji,jj,jn)
                  END DO
               END DO
            ENDIF
            !
            DO jj = 2, jpjm1                 ! Second derivative (divergence) time the eddy diffusivity coefficient
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zbtr = 1.0 / ( e12t(ji,jj) * fse3t_n(ji,jj,jk) )
                  zlt(ji,jj) = fsahtt(ji,jj,jk) * zbtr * (   ztu(ji,jj,jk) - ztu(ji-1,jj,jk)   &
                     &                                     + ztv(ji,jj,jk) - ztv(ji,jj-1,jk)   )
               END DO
            END DO
            CALL lbc_lnk( zlt, 'T', 1. )     ! Lateral boundary conditions (unchanged sgn)

            !                          !==  Bilaplacian  ==!
            !
            DO jj = 1, jpjm1                 ! third derivative (gradient)
               DO ji = 1, fs_jpim1   ! vector opt.
                  ztu(ji,jj,jk) = zeeu(ji,jj) * ( zlt(ji+1,jj  ) - zlt(ji,jj) )
                  ztv(ji,jj,jk) = zeev(ji,jj) * ( zlt(ji  ,jj+1) - zlt(ji,jj) )
               END DO
            END DO
            DO jj = 2, jpjm1                 ! fourth derivative (divergence) and add to the general tracer trend
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  ! horizontal diffusive trends
                  zbtr = 1.0 / ( e12t(ji,jj) * fse3t_n(ji,jj,jk) )
                  ztra = zbtr * (  ztu(ji,jj,jk) - ztu(ji-1,jj,jk) + ztv(ji,jj,jk) - ztv(ji,jj-1,jk)  )
                  ! add it to the general tracer trends
                  pta(ji,jj,jk,jn) = pta(ji,jj,jk,jn) + ztra
               END DO
            END DO
            !                                             
         END DO                                           ! Horizontal slab
         !                                                
         ! "zonal" mean lateral diffusive heat and salt transport
         IF( cdtype == 'TRA' .AND. ln_diaptr ) THEN  
           IF( jn == jp_tem )  htr_ldf(:) = ptr_sj( ztv(:,:,:) )
           IF( jn == jp_sal )  str_ldf(:) = ptr_sj( ztv(:,:,:) )
         ENDIF
         !
         IF( PRESENT( ptru ) ) THEN
            ptru(:,:,:,jn) = ztu(:,:,:)
            ptrv(:,:,:,jn) = ztv(:,:,:)
            ptrw(:,:,:,jn) = 0.
         ENDIF
         !                                                ! ===========
      END DO                                              ! tracer loop
      !                                                   ! ===========
      IF( nn_timing == 1 )  CALL timing_stop( 'tra_ldf_bilap')
      !
      CALL wrk_dealloc( jpi, jpj, zeeu, zeev, zlt ) 
      !
   END SUBROUTINE tra_ldf_bilap

   !!==============================================================================
END MODULE traldf_bilap
