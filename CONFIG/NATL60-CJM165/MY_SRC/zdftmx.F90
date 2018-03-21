MODULE zdftmx
   !!========================================================================
   !!                       ***  MODULE  zdftmx  ***
   !! Ocean physics: vertical tidal mixing coefficient
   !!========================================================================
   !! History :  1.0  !  2004-04  (L. Bessieres, G. Madec)  Original code
   !!             -   !  2006-08  (A. Koch-Larrouy) Indonesian strait
   !!            3.3  !  2010-10  (C. Ethe, G. Madec) reorganisation of initialisation phase
   !!----------------------------------------------------------------------
#if defined key_zdftmx   ||   defined key_esopa
   !!----------------------------------------------------------------------
   !!   'key_zdftmx'                                  Tidal vertical mixing
   !!----------------------------------------------------------------------
   !!   zdf_tmx       : global     momentum & tracer Kz with tidal induced Kz
   !!   tmx_itf       : Indonesian momentum & tracer Kz with tidal induced Kz 
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers variables
   USE dom_oce        ! ocean space and time domain variables
   USE zdf_oce        ! ocean vertical physics variables
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE eosbn2         ! ocean equation of state
   USE phycst         ! physical constants
   USE prtctl         ! Print control
   USE in_out_manager ! I/O manager
   USE iom            ! I/O Manager
   USE lib_mpp        ! MPP library
   USE wrk_nemo       ! work arrays
   USE timing         ! Timing
   USE lib_fortran    ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)  
   USE fldread        ! for file management
   IMPLICIT NONE
   PRIVATE

   PUBLIC   zdf_tmx         ! called in step module 
   PUBLIC   zdf_tmx_init    ! called in opa module 
   PUBLIC   zdf_tmx_alloc   ! called in nemogcm module

   LOGICAL, PUBLIC, PARAMETER ::   lk_zdftmx = .TRUE.    !: tidal mixing flag

   !                       !!* Namelist  namzdf_tmx : tidal mixing *
   REAL(wp) ::  rn_htmx     ! vertical decay scale for turbulence (meters)
   REAL(wp) ::  rn_n2min    ! threshold of the Brunt-Vaisala frequency (s-1)
   REAL(wp) ::  rn_tfe      ! tidal dissipation efficiency (St Laurent et al. 2002)
   REAL(wp) ::  rn_me       ! mixing efficiency (Osborn 1980)
   LOGICAL  ::  ln_tmx_itf  ! Indonesian Through Flow (ITF): Koch-Larrouy et al. (2007) parameterization
   REAL(wp) ::  rn_tfe_itf  ! ITF tidal dissipation efficiency (St Laurent et al. 2002)

   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   en_tmx     ! energy available for tidal mixing (W/m2)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   mask_itf   ! mask to use over Indonesian area
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   az_tmx     ! coefficient used to evaluate the tidal induced Kz

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO Consortium (2011)
   !! $Id: zdftmx.F90 5130 2015-03-05 19:59:13Z smasson $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION zdf_tmx_alloc()
      !!----------------------------------------------------------------------
      !!                ***  FUNCTION zdf_tmx_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE(en_tmx(jpi,jpj), mask_itf(jpi,jpj), az_tmx(jpi,jpj,jpk), STAT=zdf_tmx_alloc )
      !
      IF( lk_mpp             )   CALL mpp_sum ( zdf_tmx_alloc )
      IF( zdf_tmx_alloc /= 0 )   CALL ctl_warn('zdf_tmx_alloc: failed to allocate arrays')
   END FUNCTION zdf_tmx_alloc


   SUBROUTINE zdf_tmx( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE zdf_tmx  ***
      !!                   
      !! ** Purpose :   add to the vertical mixing coefficients the effect of
      !!              tidal mixing (Simmons et al 2004).
      !!
      !! ** Method  : - tidal-induced vertical mixing is given by:
      !!                  Kz_tides = az_tmx / max( rn_n2min, N^2 )
      !!              where az_tmx is a coefficient that specified the 3D space 
      !!              distribution of the faction of tidal energy taht is used
      !!              for mixing. Its expression is set in zdf_tmx_init routine,
      !!              following Simmons et al. 2004.
      !!                NB: a specific bounding procedure is performed on av_tide
      !!              so that the input tidal energy is actually almost used. The
      !!              basic maximum value is 60 cm2/s, but values of 300 cm2/s 
      !!              can be reached in area where bottom stratification is too 
      !!              weak.
      !!
      !!              - update av_tide in the Indonesian Through Flow area
      !!              following Koch-Larrouy et al. (2007) parameterisation
      !!              (see tmx_itf routine).
      !!
      !!              - update the model vertical eddy viscosity and diffusivity: 
      !!                     avt  = avt  +    av_tides
      !!                     avm  = avm  +    av_tides
      !!                     avmu = avmu + mi(av_tides)
      !!                     avmv = avmv + mj(av_tides)
      !!
      !! ** Action  :   avt, avm, avmu, avmv   increased by tidal mixing
      !!
      !! References : Simmons et al. 2004, Ocean Modelling, 6, 3-4, 245-263.
      !!              Koch-Larrouy et al. 2007, GRL.
      !!----------------------------------------------------------------------
      USE oce, zav_tide  =>   ua    ! use ua as workspace
      !!
      INTEGER, INTENT(in) ::   kt   ! ocean time-step 
      !!
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   ztpc         ! scalar workspace
      REAL(wp), POINTER, DIMENSION(:,:) ::   zkz
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('zdf_tmx')
      !
      CALL wrk_alloc( jpi,jpj, zkz )

      !                          ! ----------------------- !
      !                          !  Standard tidal mixing  !  (compute zav_tide)
      !                          ! ----------------------- !
      !                             !* First estimation (with n2 bound by rn_n2min) bounded by 60 cm2/s
      zav_tide(:,:,:) = MIN(  60.e-4, az_tmx(:,:,:) / MAX( rn_n2min, rn2(:,:,:) )  )

      zkz(:,:) = 0.e0               !* Associated potential energy consummed over the whole water column
      DO jk = 2, jpkm1
         zkz(:,:) = zkz(:,:) + fse3w(:,:,jk) * MAX( 0.e0, rn2(:,:,jk) ) * rau0 * zav_tide(:,:,jk) * wmask(:,:,jk)
      END DO

      DO jj = 1, jpj                !* Here zkz should be equal to en_tmx ==> multiply by en_tmx/zkz to recover en_tmx
         DO ji = 1, jpi
            IF( zkz(ji,jj) /= 0.e0 )   zkz(ji,jj) = en_tmx(ji,jj) / zkz(ji,jj)
         END DO
      END DO

      DO jk = 2, jpkm1     !* Mutiply by zkz to recover en_tmx, BUT bound by 30/6 ==> zav_tide bound by 300 cm2/s
         DO jj = 1, jpj                !* Here zkz should be equal to en_tmx ==> multiply by en_tmx/zkz to recover en_tmx
            DO ji = 1, jpi
               zav_tide(ji,jj,jk) = zav_tide(ji,jj,jk) * MIN( zkz(ji,jj), 30./6. ) * wmask(ji,jj,jk)  !kz max = 300 cm2/s
            END DO
         END DO
      END DO

      IF( kt == nit000 ) THEN       !* check at first time-step: diagnose the energy consumed by zav_tide
         ztpc = 0.e0
         DO jk= 1, jpk
            DO jj= 1, jpj
               DO ji= 1, jpi
                  ztpc = ztpc + fse3w(ji,jj,jk) * e1t(ji,jj) * e2t(ji,jj)   &
                     &         * MAX( 0.e0, rn2(ji,jj,jk) ) * zav_tide(ji,jj,jk) * tmask(ji,jj,jk) * tmask_i(ji,jj)
               END DO
            END DO
         END DO
         ztpc= rau0 / ( rn_tfe * rn_me ) * ztpc
         IF(lwp) WRITE(numout,*) 
         IF(lwp) WRITE(numout,*) '          N Total power consumption by av_tide    : ztpc = ', ztpc * 1.e-12 ,'TW'
      ENDIF
       
      !                          ! ----------------------- !
      !                          !    ITF  tidal mixing    !  (update zav_tide)
      !                          ! ----------------------- !
      IF( ln_tmx_itf )   CALL tmx_itf( kt, zav_tide )

      !                          ! ----------------------- !
      !                          !   Update  mixing coefs  !                          
      !                          ! ----------------------- !
      DO jk = 2, jpkm1              !* update momentum & tracer diffusivity with tidal mixing
         DO jj = 1, jpj                !* Here zkz should be equal to en_tmx ==> multiply by en_tmx/zkz to recover en_tmx
            DO ji = 1, jpi
               avt(ji,jj,jk) = avt(ji,jj,jk) + zav_tide(ji,jj,jk) * wmask(ji,jj,jk)
               avm(ji,jj,jk) = avm(ji,jj,jk) + zav_tide(ji,jj,jk) * wmask(ji,jj,jk)
            END DO
         END DO
      END DO
      
      DO jk = 2, jpkm1              !* update momentum & tracer diffusivity with tidal mixing
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1  ! vector opt.
               avmu(ji,jj,jk) = avmu(ji,jj,jk) + 0.5 * ( zav_tide(ji,jj,jk) + zav_tide(ji+1,jj  ,jk) ) * wumask(ji,jj,jk)
               avmv(ji,jj,jk) = avmv(ji,jj,jk) + 0.5 * ( zav_tide(ji,jj,jk) + zav_tide(ji  ,jj+1,jk) ) * wvmask(ji,jj,jk)
            END DO
         END DO
      END DO
      CALL lbc_lnk( avmu, 'U', 1. )   ;   CALL lbc_lnk( avmv, 'V', 1. )      ! lateral boundary condition

      !                             !* output tidal mixing coefficient
      CALL iom_put( "av_tide", zav_tide )

      IF(ln_ctl)   CALL prt_ctl(tab3d_1=zav_tide , clinfo1=' tmx - av_tide: ', tab3d_2=avt, clinfo2=' avt: ', ovlap=1, kdim=jpk)
      !
      CALL wrk_dealloc( jpi,jpj, zkz )
      !
      IF( nn_timing == 1 )  CALL timing_stop('zdf_tmx')
      !
   END SUBROUTINE zdf_tmx


   SUBROUTINE tmx_itf( kt, pav )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tmx_itf  ***
      !!                   
      !! ** Purpose :   modify the vertical eddy diffusivity coefficients 
      !!              (pav) in the Indonesian Through Flow area (ITF).
      !!
      !! ** Method  : - Following Koch-Larrouy et al. (2007), in the ITF defined
      !!                by msk_itf (read in a file, see tmx_init), the tidal
      !!                mixing coefficient is computed with :
      !!                  * q=1 (i.e. all the tidal energy remains trapped in
      !!                         the area and thus is used for mixing)
      !!                  * the vertical distribution of the tifal energy is a
      !!                    proportional to N above the thermocline (d(N^2)/dz > 0)
      !!                    and to N^2 below the thermocline (d(N^2)/dz < 0)
      !!
      !! ** Action  :   av_tide   updated in the ITF area (msk_itf)
      !!
      !! References :  Koch-Larrouy et al. 2007, GRL 
      !!----------------------------------------------------------------------
      INTEGER , INTENT(in   )                         ::   kt   ! ocean time-step
      REAL(wp), INTENT(inout), DIMENSION(jpi,jpj,jpk) ::   pav  ! Tidal mixing coef.
      !! 
      INTEGER  ::   ji, jj, jk    ! dummy loop indices
      REAL(wp) ::   zcoef, ztpc   ! temporary scalar
      REAL(wp), DIMENSION(:,:)  , POINTER ::   zkz                        ! 2D workspace
      REAL(wp), DIMENSION(:,:)  , POINTER ::   zsum1 , zsum2 , zsum       !  -      -
      REAL(wp), DIMENSION(:,:,:), POINTER ::   zempba_3d_1, zempba_3d_2   ! 3D workspace
      REAL(wp), DIMENSION(:,:,:), POINTER ::   zempba_3d  , zdn2dz        !  -      -
      REAL(wp), DIMENSION(:,:,:), POINTER ::   zavt_itf                   !  -      -
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('tmx_itf')
      !
      CALL wrk_alloc( jpi,jpj, zkz, zsum1 , zsum2 , zsum )
      CALL wrk_alloc( jpi,jpj,jpk, zempba_3d_1, zempba_3d_2, zempba_3d, zdn2dz, zavt_itf )

      !                             ! compute the form function using N2 at each time step
      zempba_3d_1(:,:,jpk) = 0.e0
      zempba_3d_2(:,:,jpk) = 0.e0
      DO jk = 1, jpkm1             
         zdn2dz     (:,:,jk) = rn2(:,:,jk) - rn2(:,:,jk+1)           ! Vertical profile of dN2/dz
!CDIR NOVERRCHK
         zempba_3d_1(:,:,jk) = SQRT(  MAX( 0.e0, rn2(:,:,jk) )  )    !    -        -    of N
         zempba_3d_2(:,:,jk) =        MAX( 0.e0, rn2(:,:,jk) )       !    -        -    of N^2
      END DO
      !
      zsum (:,:) = 0.e0
      zsum1(:,:) = 0.e0
      zsum2(:,:) = 0.e0
      DO jk= 2, jpk
         zsum1(:,:) = zsum1(:,:) + zempba_3d_1(:,:,jk) * fse3w(:,:,jk) * tmask(:,:,jk) * tmask(:,:,jk-1)
         zsum2(:,:) = zsum2(:,:) + zempba_3d_2(:,:,jk) * fse3w(:,:,jk) * tmask(:,:,jk) * tmask(:,:,jk-1)               
      END DO
      DO jj = 1, jpj
         DO ji = 1, jpi
            IF( zsum1(ji,jj) /= 0.e0 )   zsum1(ji,jj) = 1.e0 / zsum1(ji,jj)
            IF( zsum2(ji,jj) /= 0.e0 )   zsum2(ji,jj) = 1.e0 / zsum2(ji,jj)                
         END DO
      END DO

      DO jk= 1, jpk
         DO jj = 1, jpj
            DO ji = 1, jpi
               zcoef = 0.5 - SIGN( 0.5, zdn2dz(ji,jj,jk) )       ! =0 if dN2/dz > 0, =1 otherwise 
               ztpc  = zempba_3d_1(ji,jj,jk) * zsum1(ji,jj) *        zcoef     &
                  &  + zempba_3d_2(ji,jj,jk) * zsum2(ji,jj) * ( 1. - zcoef )
               !
               zempba_3d(ji,jj,jk) =               ztpc 
               zsum     (ji,jj)    = zsum(ji,jj) + ztpc * fse3w(ji,jj,jk)
            END DO
         END DO
       END DO
       DO jj = 1, jpj
          DO ji = 1, jpi
             IF( zsum(ji,jj) > 0.e0 )   zsum(ji,jj) = 1.e0 / zsum(ji,jj)                
          END DO
       END DO

      !                             ! first estimation bounded by 10 cm2/s (with n2 bounded by rn_n2min) 
      zcoef = rn_tfe_itf / ( rn_tfe * rau0 )
      DO jk = 1, jpk
         zavt_itf(:,:,jk) = MIN(  10.e-4, zcoef * en_tmx(:,:) * zsum(:,:) * zempba_3d(:,:,jk)   &
            &                                      / MAX( rn_n2min, rn2(:,:,jk) ) * tmask(:,:,jk)  )
      END DO           

      zkz(:,:) = 0.e0               ! Associated potential energy consummed over the whole water column
      DO jk = 2, jpkm1
         zkz(:,:) = zkz(:,:) + fse3w(:,:,jk) * MAX( 0.e0, rn2(:,:,jk) ) * rau0 * zavt_itf(:,:,jk) * tmask(:,:,jk) * tmask(:,:,jk-1)
      END DO

      DO jj = 1, jpj                ! Here zkz should be equal to en_tmx ==> multiply by en_tmx/zkz to recover en_tmx
         DO ji = 1, jpi
            IF( zkz(ji,jj) /= 0.e0 )   zkz(ji,jj) = en_tmx(ji,jj) * rn_tfe_itf / rn_tfe / zkz(ji,jj)
         END DO
      END DO

      DO jk = 2, jpkm1              ! Mutiply by zkz to recover en_tmx, BUT bound by 30/6 ==> zavt_itf bound by 300 cm2/s
         zavt_itf(:,:,jk) = zavt_itf(:,:,jk) * MIN( zkz(:,:), 120./10. ) * tmask(:,:,jk) * tmask(:,:,jk-1)   ! kz max = 120 cm2/s
      END DO

      IF( kt == nit000 ) THEN       ! diagnose the nergy consumed by zavt_itf
         ztpc = 0.e0
         DO jk= 1, jpk
            DO jj= 1, jpj
               DO ji= 1, jpi
                  ztpc = ztpc + e1t(ji,jj) * e2t(ji,jj) * fse3w(ji,jj,jk) * MAX( 0.e0, rn2(ji,jj,jk) )   &
                     &                     * zavt_itf(ji,jj,jk) * tmask(ji,jj,jk) * tmask_i(ji,jj)
               END DO
            END DO
         END DO
         ztpc= rau0 * ztpc / ( rn_me * rn_tfe_itf )
         IF(lwp) WRITE(numout,*) '          N Total power consumption by zavt_itf: ztpc = ', ztpc * 1.e-12 ,'TW'
      ENDIF

      !                             ! Update pav with the ITF mixing coefficient
      DO jk = 2, jpkm1
         pav(:,:,jk) = pav     (:,:,jk) * ( 1.e0 - mask_itf(:,:) )   &
            &        + zavt_itf(:,:,jk) *          mask_itf(:,:) 
      END DO
      !
      CALL wrk_dealloc( jpi,jpj, zkz, zsum1 , zsum2 , zsum )
      CALL wrk_dealloc( jpi,jpj,jpk, zempba_3d_1, zempba_3d_2, zempba_3d, zdn2dz, zavt_itf )
      !
      IF( nn_timing == 1 )  CALL timing_stop('tmx_itf')
      !
   END SUBROUTINE tmx_itf


   SUBROUTINE zdf_tmx_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE zdf_tmx_init  ***
      !!                     
      !! ** Purpose :   Initialization of the vertical tidal mixing, Reading
      !!              of M2 and K1 tidal energy in nc files
      !!
      !! ** Method  : - Read the namtmx namelist and check the parameters
      !!
      !!              - Read the input data in NetCDF files :
      !!              M2 and K1 tidal energy. The total tidal energy, en_tmx, 
      !!              is the sum of M2, K1 and S2 energy where S2 is assumed 
      !!              to be: S2=(1/2)^2 * M2
      !!              mask_itf, a mask array that determine where substituing 
      !!              the standard Simmons et al. (2005) formulation with the
      !!              one of Koch_Larrouy et al. (2007).
      !!
      !!              - Compute az_tmx, a 3D coefficient that allows to compute
      !!             the standard tidal-induced vertical mixing as follows:
      !!                  Kz_tides = az_tmx / max( rn_n2min, N^2 )
      !!             with az_tmx a bottom intensified coefficient is given by:
      !!                 az_tmx(z) = en_tmx / ( rau0 * rn_htmx ) * EXP( -(H-z)/rn_htmx )
      !!                                                  / ( 1. - EXP( - H   /rn_htmx ) ) 
      !!             where rn_htmx the characteristic length scale of the bottom 
      !!             intensification, en_tmx the tidal energy, and H the ocean depth
      !!
      !! ** input   :   - Namlist namtmx
      !!                - NetCDF file : M2_ORCA2.nc, K1_ORCA2.nc, and mask_itf.nc
      !!
      !! ** Action  : - Increase by 1 the nstop flag is setting problem encounter
      !!              - defined az_tmx used to compute tidal-induced mixing
      !!
      !! References : Simmons et al. 2004, Ocean Modelling, 6, 3-4, 245-263.
      !!              Koch-Larrouy et al. 2007, GRL.
      !!----------------------------------------------------------------------
      USE oce     ,         zav_tide =>  ua         ! ua used as workspace
      !!
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      INTEGER  ::   inum         ! local integer
      INTEGER  ::   ios
      REAL(wp) ::   ztpc, ze_z   ! local scalars
      REAL(wp), DIMENSION(:,:)  , POINTER ::  zem2, zek1 ,zes2   ! read M2 and K1 tidal energy
      REAL(wp), DIMENSION(:,:)  , POINTER ::  zkz                ! total M2, K1 and S2 tidal energy
      REAL(wp), DIMENSION(:,:)  , POINTER ::  zfact              ! used for vertical structure function
      REAL(wp), DIMENSION(:,:)  , POINTER ::  zhdep              ! Ocean depth 
      REAL(wp), DIMENSION(:,:,:), POINTER ::  zpc                ! power consumption
      TYPE(FLD_N) :: sn_mskitf, sn_m2, sn_s2, sn_k1
      !!
      NAMELIST/namzdf_tmx/ rn_htmx, rn_n2min, rn_tfe, rn_me, ln_tmx_itf, rn_tfe_itf,     &
         &    sn_mskitf, sn_m2, sn_s2, sn_k1
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('zdf_tmx_init')
      !
      CALL wrk_alloc( jpi,jpj, zem2, zes2, zek1, zkz, zfact, zhdep )
      CALL wrk_alloc( jpi,jpj,jpk, zpc )
      
      REWIND( numnam_ref )              ! Namelist namzdf_tmx in reference namelist : Tidal Mixing
      READ  ( numnam_ref, namzdf_tmx, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namzdf_tmx in reference namelist', lwp )

      REWIND( numnam_cfg )              ! Namelist namzdf_tmx in configuration namelist : Tidal Mixing
      READ  ( numnam_cfg, namzdf_tmx, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namzdf_tmx in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, namzdf_tmx )

      IF(lwp) THEN                   ! Control print
         WRITE(numout,*)
         WRITE(numout,*) 'zdf_tmx_init : tidal mixing'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namzdf_tmx : set tidal mixing parameters'
         WRITE(numout,*) '      Vertical decay scale for turbulence   = ', rn_htmx 
         WRITE(numout,*) '      Brunt-Vaisala frequency threshold     = ', rn_n2min
         WRITE(numout,*) '      Tidal dissipation efficiency          = ', rn_tfe
         WRITE(numout,*) '      Mixing efficiency                     = ', rn_me
         WRITE(numout,*) '      ITF specific parameterisation         = ', ln_tmx_itf
         WRITE(numout,*) '      ITF tidal dissipation efficiency      = ', rn_tfe_itf
      ENDIF

      !                              ! allocate tmx arrays
      IF( zdf_tmx_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'zdf_tmx_init : unable to allocate tmx arrays' )

      IF( ln_tmx_itf ) THEN          ! read the Indonesian Through Flow mask
         CALL iom_open  (sn_mskitf%clname, inum)
         CALL iom_get   (inum, jpdom_data, sn_mskitf%clvar, mask_itf,1) ! 
         CALL iom_close (inum)
      ENDIF

      ! read M2 tidal energy flux : W/m2  ( zem2 < 0 )
      CALL iom_open  (sn_m2%clname, inum)
      CALL iom_get   (inum, jpdom_data, sn_m2%clvar,  zem2,1) ! 
      WHERE ( zem2 > 0.d0 ) zem2 = 0.d0
      CALL iom_close (inum)

      ! read S2 tidal energy flux : W/m2  ( zes2 < 0 )
      CALL iom_open  (sn_s2%clname, inum)
      CALL iom_get   (inum, jpdom_data, sn_s2%clvar,  zes2,1) ! 
      WHERE ( zes2 > 0.d0 ) zes2 = 0.d0
      CALL iom_close (inum)

      ! read K1 tidal energy flux : W/m2  ( zek1 < 0 )
      CALL iom_open  (sn_k1%clname, inum)
      CALL iom_get   (inum, jpdom_data, sn_k1%clvar,  zek1,1) ! 
      WHERE ( zek1 > 0.d0 ) zek1 = 0.d0
      CALL iom_close (inum)
 
      ! Total tidal energy ( M2, S2 and K1  with S2=(1/2)^2 * M2 )
      ! only the energy available for mixing is taken into account,
      ! (mixing efficiency tidal dissipation efficiency)
      !  zes2(:,:) = zem2(:,:) * 0.25   !  if no specific field for S2
      en_tmx(:,:) = - rn_tfe * rn_me * ( zem2(:,:) + zes2(:,:) + zek1(:,:) ) * ssmask(:,:)

!============
!TG: Bug for VVL? Should this section be moved out of _init and be updated at every timestep?
      ! Vertical structure (az_tmx)
      DO jj = 1, jpj                ! part independent of the level
         DO ji = 1, jpi
            zhdep(ji,jj) = gdepw_0(ji,jj,mbkt(ji,jj)+1)       ! depth of the ocean
            zfact(ji,jj) = rau0 * rn_htmx * ( 1. - EXP( -zhdep(ji,jj) / rn_htmx ) )
            IF( zfact(ji,jj) /= 0 )   zfact(ji,jj) = en_tmx(ji,jj) / zfact(ji,jj)
         END DO
      END DO
      DO jk= 1, jpk                 ! complete with the level-dependent part
         DO jj = 1, jpj
            DO ji = 1, jpi
               az_tmx(ji,jj,jk) = zfact(ji,jj) * EXP( -( zhdep(ji,jj)-gdepw_0(ji,jj,jk) ) / rn_htmx ) * tmask(ji,jj,jk)
            END DO
         END DO
      END DO
!===========

      IF( nprint == 1 .AND. lwp ) THEN
         ! Control print
         ! Total power consumption due to vertical mixing
         ! zpc = rau0 * 1/rn_me * rn2 * zav_tide
         zav_tide(:,:,:) = 0.e0
         DO jk = 2, jpkm1
            zav_tide(:,:,jk) = az_tmx(:,:,jk) / MAX( rn_n2min, rn2(:,:,jk) )
         END DO

         ztpc = 0.e0
         zpc(:,:,:) = MAX(rn_n2min,rn2(:,:,:)) * zav_tide(:,:,:)
         DO jk= 2, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  ztpc = ztpc + fse3w(ji,jj,jk) * e1t(ji,jj) * e2t(ji,jj) * zpc(ji,jj,jk) * wmask(ji,jj,jk) * tmask_i(ji,jj)
               END DO
            END DO
         END DO
         ztpc= rau0 * 1/(rn_tfe * rn_me) * ztpc

         WRITE(numout,*) 
         WRITE(numout,*) '          Total power consumption of the tidally driven part of Kz : ztpc = ', ztpc * 1.e-12 ,'TW'


         ! control print 2
         zav_tide(:,:,:) = MIN( zav_tide(:,:,:), 60.e-4 )   
         zkz(:,:) = 0.e0
         DO jk = 2, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zkz(ji,jj) = zkz(ji,jj) + fse3w(ji,jj,jk) * MAX(0.e0, rn2(ji,jj,jk)) * rau0 * zav_tide(ji,jj,jk) * wmask(ji,jj,jk)
               END DO
            END DO
         END DO
         ! Here zkz should be equal to en_tmx ==> multiply by en_tmx/zkz
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( zkz(ji,jj) /= 0.e0 )   THEN
                   zkz(ji,jj) = en_tmx(ji,jj) / zkz(ji,jj)
               ENDIF
            END DO
         END DO
         ztpc = 1.e50
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( zkz(ji,jj) /= 0.e0 )   THEN
                   ztpc = Min( zkz(ji,jj), ztpc)
               ENDIF
            END DO
         END DO
         WRITE(numout,*) '          Min de zkz ', ztpc, ' Max = ', maxval(zkz(:,:) )

         DO jk = 2, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zav_tide(ji,jj,jk) = zav_tide(ji,jj,jk) * MIN( zkz(ji,jj), 30./6. ) * wmask(ji,jj,jk)  !kz max = 300 cm2/s
               END DO
            END DO
         END DO
         ztpc = 0.e0
         zpc(:,:,:) = Max(0.e0,rn2(:,:,:)) * zav_tide(:,:,:)
         DO jk= 1, jpk
            DO jj = 1, jpj
               DO ji = 1, jpi
                  ztpc = ztpc + fse3w(ji,jj,jk) * e1t(ji,jj) * e2t(ji,jj) * zpc(ji,jj,jk) * wmask(ji,jj,jk) * tmask_i(ji,jj)
               END DO
            END DO
         END DO
         ztpc= rau0 * 1/(rn_tfe * rn_me) * ztpc
         WRITE(numout,*) '          2 Total power consumption of the tidally driven part of Kz : ztpc = ', ztpc * 1.e-12 ,'TW'

         DO jk = 1, jpk
            ze_z =                  SUM( e1t(:,:) * e2t(:,:) * zav_tide(:,:,jk)     * tmask_i(:,:) )   &
               &     / MAX( 1.e-20, SUM( e1t(:,:) * e2t(:,:) * wmask (:,:,jk) * tmask_i(:,:) ) )
            ztpc = 1.E50
            DO jj = 1, jpj
               DO ji = 1, jpi
                  IF( zav_tide(ji,jj,jk) /= 0.e0 )   ztpc =Min( ztpc, zav_tide(ji,jj,jk) )
               END DO
            END DO
            WRITE(numout,*) '            N2 min - jk= ', jk,'   ', ze_z * 1.e4,' cm2/s min= ',ztpc*1.e4,   &
               &       'max= ', MAXVAL(zav_tide(:,:,jk) )*1.e4, ' cm2/s'
         END DO

         WRITE(numout,*) '          e_tide : ', SUM( e1t*e2t*en_tmx ) / ( rn_tfe * rn_me ) * 1.e-12, 'TW'
         WRITE(numout,*) 
         WRITE(numout,*) '          Initial profile of tidal vertical mixing'
         DO jk = 1, jpk
            DO jj = 1,jpj
               DO ji = 1,jpi
                  zkz(ji,jj) = az_tmx(ji,jj,jk) /MAX( rn_n2min, rn2(ji,jj,jk) )
               END DO
            END DO
            ze_z =                  SUM( e1t(:,:) * e2t(:,:) * zkz(:,:)     * tmask_i(:,:) )   &
               &     / MAX( 1.e-20, SUM( e1t(:,:) * e2t(:,:) * wmask (:,:,jk) * tmask_i(:,:) ) )
            WRITE(numout,*) '                jk= ', jk,'   ', ze_z * 1.e4,' cm2/s'
         END DO
         DO jk = 1, jpk
            zkz(:,:) = az_tmx(:,:,jk) /rn_n2min
            ze_z =                  SUM( e1t(:,:) * e2t(:,:) * zkz(:,:)     * tmask_i(:,:) )   &
               &     / MAX( 1.e-20, SUM( e1t(:,:) * e2t(:,:) * wmask (:,:,jk) * tmask_i(:,:) ) )
            WRITE(numout,*) 
            WRITE(numout,*) '          N2 min - jk= ', jk,'   ', ze_z * 1.e4,' cm2/s min= ',MINVAL(zkz)*1.e4,   &
               &       'max= ', MAXVAL(zkz)*1.e4, ' cm2/s'
         END DO
         !
      ENDIF
      !
      CALL wrk_dealloc( jpi,jpj, zem2, zek1, zes2, zkz, zfact, zhdep )
      CALL wrk_dealloc( jpi,jpj,jpk, zpc )
      !
      IF( nn_timing == 1 )  CALL timing_stop('zdf_tmx_init')
      !
   END SUBROUTINE zdf_tmx_init

#elif defined key_zdftmx_new
   !!----------------------------------------------------------------------
   !!   'key_zdftmx_new'               Internal wave-driven vertical mixing
   !!----------------------------------------------------------------------
   !!   zdf_tmx       : global     momentum & tracer Kz with wave induced Kz
   !!   zdf_tmx_init  : global     momentum & tracer Kz with wave induced Kz
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers variables
   USE dom_oce        ! ocean space and time domain variables
   USE zdf_oce        ! ocean vertical physics variables
   USE zdfddm         ! ocean vertical physics: double diffusive mixing
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE eosbn2         ! ocean equation of state
   USE phycst         ! physical constants
   USE prtctl         ! Print control
   USE in_out_manager ! I/O manager
   USE iom            ! I/O Manager
   USE lib_mpp        ! MPP library
   USE wrk_nemo       ! work arrays
   USE timing         ! Timing
   USE lib_fortran    ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)  

   IMPLICIT NONE
   PRIVATE

   PUBLIC   zdf_tmx         ! called in step module 
   PUBLIC   zdf_tmx_init    ! called in nemogcm module 
   PUBLIC   zdf_tmx_alloc   ! called in nemogcm module

   LOGICAL, PUBLIC, PARAMETER ::   lk_zdftmx = .TRUE.    !: wave-driven mixing flag

   !                       !!* Namelist  namzdf_tmx : internal wave-driven mixing *
   INTEGER  ::  nn_zpyc     ! pycnocline-intensified mixing energy proportional to N (=1) or N^2 (=2)
   LOGICAL  ::  ln_mevar    ! variable (=T) or constant (=F) mixing efficiency
   LOGICAL  ::  ln_tsdiff   ! account for differential T/S wave-driven mixing (=T) or not (=F)

   REAL(wp) ::  r1_6 = 1._wp / 6._wp

   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   ebot_tmx     ! power available from high-mode wave breaking (W/m2)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   epyc_tmx     ! power available from low-mode, pycnocline-intensified wave breaking (W/m2)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   ecri_tmx     ! power available from low-mode, critical slope wave breaking (W/m2)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   hbot_tmx     ! WKB decay scale for high-mode energy dissipation (m)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   hcri_tmx     ! decay scale for low-mode critical slope dissipation (m)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   emix_tmx     ! local energy density available for mixing (W/kg)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   bflx_tmx     ! buoyancy flux Kz * N^2 (W/kg)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   pcmap_tmx    ! vertically integrated buoyancy flux (W/m2)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   zav_ratio    ! S/T diffusivity ratio (only for ln_tsdiff=T)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   zav_wave     ! Internal wave-induced diffusivity

   !! * Substitutions
#  include "zdfddm_substitute.h90"
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO Consortium (2016)
   !! $Id: zdftmx.F90 6314 2016-02-15 12:04:56Z cetlod $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION zdf_tmx_alloc()
      !!----------------------------------------------------------------------
      !!                ***  FUNCTION zdf_tmx_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE(     ebot_tmx(jpi,jpj),  epyc_tmx(jpi,jpj),  ecri_tmx(jpi,jpj)    ,   &
      &             hbot_tmx(jpi,jpj),  hcri_tmx(jpi,jpj),  emix_tmx(jpi,jpj,jpk),   &
      &         bflx_tmx(jpi,jpj,jpk), pcmap_tmx(jpi,jpj), zav_ratio(jpi,jpj,jpk),   & 
      &         zav_wave(jpi,jpj,jpk), STAT=zdf_tmx_alloc     )
      !
      IF( lk_mpp             )   CALL mpp_sum ( zdf_tmx_alloc )
      IF( zdf_tmx_alloc /= 0 )   CALL ctl_warn('zdf_tmx_alloc: failed to allocate arrays')
   END FUNCTION zdf_tmx_alloc


   SUBROUTINE zdf_tmx( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE zdf_tmx  ***
      !!                   
      !! ** Purpose :   add to the vertical mixing coefficients the effect of
      !!              breaking internal waves.
      !!
      !! ** Method  : - internal wave-driven vertical mixing is given by:
      !!                  Kz_wave = min(  100 cm2/s, f(  Reb = emix_tmx /( Nu * N^2 )  )
      !!              where emix_tmx is the 3D space distribution of the wave-breaking 
      !!              energy and Nu the molecular kinematic viscosity.
      !!              The function f(Reb) is linear (constant mixing efficiency)
      !!              if the namelist parameter ln_mevar = F and nonlinear if ln_mevar = T.
      !!
      !!              - Compute emix_tmx, the 3D power density that allows to compute
      !!              Reb and therefrom the wave-induced vertical diffusivity.
      !!              This is divided into three components:
      !!                 1. Bottom-intensified low-mode dissipation at critical slopes
      !!                     emix_tmx(z) = ( ecri_tmx / rau0 ) * EXP( -(H-z)/hcri_tmx )
      !!                                   / ( 1. - EXP( - H/hcri_tmx ) ) * hcri_tmx
      !!              where hcri_tmx is the characteristic length scale of the bottom 
      !!              intensification, ecri_tmx a map of available power, and H the ocean depth.
      !!                 2. Pycnocline-intensified low-mode dissipation
      !!                     emix_tmx(z) = ( epyc_tmx / rau0 ) * ( sqrt(rn2(z))^nn_zpyc )
      !!                                   / SUM( sqrt(rn2(z))^nn_zpyc * e3w(z) )
      !!              where epyc_tmx is a map of available power, and nn_zpyc
      !!              is the chosen stratification-dependence of the internal wave
      !!              energy dissipation.
      !!                 3. WKB-height dependent high mode dissipation
      !!                     emix_tmx(z) = ( ebot_tmx / rau0 ) * rn2(z) * EXP(-z_wkb(z)/hbot_tmx)
      !!                                   / SUM( rn2(z) * EXP(-z_wkb(z)/hbot_tmx) * e3w(z) )
      !!              where hbot_tmx is the characteristic length scale of the WKB bottom 
      !!              intensification, ebot_tmx is a map of available power, and z_wkb is the
      !!              WKB-stretched height above bottom defined as
      !!                    z_wkb(z) = H * SUM( sqrt(rn2(z'>=z)) * e3w(z'>=z) )
      !!                                 / SUM( sqrt(rn2(z'))    * e3w(z')    )
      !!
      !!              - update the model vertical eddy viscosity and diffusivity: 
      !!                     avt  = avt  +    av_wave
      !!                     avm  = avm  +    av_wave
      !!                     avmu = avmu + mi(av_wave)
      !!                     avmv = avmv + mj(av_wave)
      !!
      !!              - if namelist parameter ln_tsdiff = T, account for differential mixing:
      !!                     avs  = avt  +    av_wave * diffusivity_ratio(Reb)
      !!
      !! ** Action  : - Define emix_tmx used to compute internal wave-induced mixing
      !!              - avt, avs, avm, avmu, avmv increased by internal wave-driven mixing    
      !!
      !! References :  de Lavergne et al. 2015, JPO; 2016, in prep.
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time-step 
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   ztpc         ! scalar workspace
      REAL(wp), DIMENSION(:,:)  , POINTER ::  zfact     ! Used for vertical structure
      REAL(wp), DIMENSION(:,:)  , POINTER ::  zhdep     ! Ocean depth
      REAL(wp), DIMENSION(:,:,:), POINTER ::  zwkb      ! WKB-stretched height above bottom
      REAL(wp), DIMENSION(:,:,:), POINTER ::  zweight   ! Weight for high mode vertical distribution
      REAL(wp), DIMENSION(:,:,:), POINTER ::  znu_t     ! Molecular kinematic viscosity (T grid)
      REAL(wp), DIMENSION(:,:,:), POINTER ::  znu_w     ! Molecular kinematic viscosity (W grid)
      REAL(wp), DIMENSION(:,:,:), POINTER ::  zReb      ! Turbulence intensity parameter
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )   CALL timing_start('zdf_tmx')
      !
      CALL wrk_alloc( jpi,jpj,       zfact, zhdep )
      CALL wrk_alloc( jpi,jpj,jpk,   zwkb, zweight, znu_t, znu_w, zReb )

      !                          ! ----------------------------- !
      !                          !  Internal wave-driven mixing  !  (compute zav_wave)
      !                          ! ----------------------------- !
      !                             
      !                        !* Critical slope mixing: distribute energy over the time-varying ocean depth,
      !                                                 using an exponential decay from the seafloor.
      DO jj = 1, jpj                ! part independent of the level
         DO ji = 1, jpi
            zhdep(ji,jj) = fsdepw(ji,jj,mbkt(ji,jj)+1)       ! depth of the ocean
            zfact(ji,jj) = rau0 * (  1._wp - EXP( -zhdep(ji,jj) / hcri_tmx(ji,jj) )  )
            IF( zfact(ji,jj) /= 0 )   zfact(ji,jj) = ecri_tmx(ji,jj) / zfact(ji,jj)
         END DO
      END DO

      DO jk = 2, jpkm1              ! complete with the level-dependent part
         emix_tmx(:,:,jk) = zfact(:,:) * (  EXP( ( fsde3w(:,:,jk  ) - zhdep(:,:) ) / hcri_tmx(:,:) )                      &
            &                             - EXP( ( fsde3w(:,:,jk-1) - zhdep(:,:) ) / hcri_tmx(:,:) )  ) * wmask(:,:,jk)   &
            &                          / ( fsde3w(:,:,jk) - fsde3w(:,:,jk-1) )
      END DO

      !                        !* Pycnocline-intensified mixing: distribute energy over the time-varying 
      !                        !* ocean depth as proportional to sqrt(rn2)^nn_zpyc

      SELECT CASE ( nn_zpyc )

      CASE ( 1 )               ! Dissipation scales as N (recommended)

         zfact(:,:) = 0._wp
         DO jk = 2, jpkm1              ! part independent of the level
            zfact(:,:) = zfact(:,:) + fse3w(:,:,jk) * SQRT(  MAX( 0._wp, rn2(:,:,jk) )  ) * wmask(:,:,jk)
         END DO

         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( zfact(ji,jj) /= 0 )   zfact(ji,jj) = epyc_tmx(ji,jj) / ( rau0 * zfact(ji,jj) )
            END DO
         END DO

         DO jk = 2, jpkm1              ! complete with the level-dependent part
            emix_tmx(:,:,jk) = emix_tmx(:,:,jk) + zfact(:,:) * SQRT(  MAX( 0._wp, rn2(:,:,jk) )  ) * wmask(:,:,jk)
         END DO

      CASE ( 2 )               ! Dissipation scales as N^2

         zfact(:,:) = 0._wp
         DO jk = 2, jpkm1              ! part independent of the level
            zfact(:,:) = zfact(:,:) + fse3w(:,:,jk) * MAX( 0._wp, rn2(:,:,jk) ) * wmask(:,:,jk)
         END DO

         DO jj= 1, jpj
            DO ji = 1, jpi
               IF( zfact(ji,jj) /= 0 )   zfact(ji,jj) = epyc_tmx(ji,jj) / ( rau0 * zfact(ji,jj) )
            END DO
         END DO

         DO jk = 2, jpkm1              ! complete with the level-dependent part
            emix_tmx(:,:,jk) = emix_tmx(:,:,jk) + zfact(:,:) * MAX( 0._wp, rn2(:,:,jk) ) * wmask(:,:,jk)
         END DO

      END SELECT

      !                        !* WKB-height dependent mixing: distribute energy over the time-varying 
      !                        !* ocean depth as proportional to rn2 * exp(-z_wkb/rn_hbot)
      
      zwkb(:,:,:) = 0._wp
      zfact(:,:) = 0._wp
      DO jk = 2, jpkm1
         zfact(:,:) = zfact(:,:) + fse3w(:,:,jk) * SQRT(  MAX( 0._wp, rn2(:,:,jk) )  ) * wmask(:,:,jk)
         zwkb(:,:,jk) = zfact(:,:)
      END DO

      DO jk = 2, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( zfact(ji,jj) /= 0 )   zwkb(ji,jj,jk) = zhdep(ji,jj) * ( zfact(ji,jj) - zwkb(ji,jj,jk) )   &
                                            &           * tmask(ji,jj,jk) / zfact(ji,jj)
            END DO
         END DO
      END DO
      zwkb(:,:,1) = zhdep(:,:) * tmask(:,:,1)

      zweight(:,:,:) = 0._wp
      DO jk = 2, jpkm1
         zweight(:,:,jk) = MAX( 0._wp, rn2(:,:,jk) ) * hbot_tmx(:,:) * wmask(:,:,jk)                    &
            &   * (  EXP( -zwkb(:,:,jk) / hbot_tmx(:,:) ) - EXP( -zwkb(:,:,jk-1) / hbot_tmx(:,:) )  )
      END DO

      zfact(:,:) = 0._wp
      DO jk = 2, jpkm1              ! part independent of the level
         zfact(:,:) = zfact(:,:) + zweight(:,:,jk)
      END DO

      DO jj = 1, jpj
         DO ji = 1, jpi
            IF( zfact(ji,jj) /= 0 )   zfact(ji,jj) = ebot_tmx(ji,jj) / ( rau0 * zfact(ji,jj) )
         END DO
      END DO

      DO jk = 2, jpkm1              ! complete with the level-dependent part
         emix_tmx(:,:,jk) = emix_tmx(:,:,jk) + zweight(:,:,jk) * zfact(:,:) * wmask(:,:,jk)   &
            &                                / ( fsde3w(:,:,jk) - fsde3w(:,:,jk-1) )
      END DO


      ! Calculate molecular kinematic viscosity
      znu_t(:,:,:) = 1.e-4_wp * (  17.91_wp - 0.53810_wp * tsn(:,:,:,jp_tem) + 0.00694_wp * tsn(:,:,:,jp_tem) * tsn(:,:,:,jp_tem)  &
         &                                  + 0.02305_wp * tsn(:,:,:,jp_sal)  ) * tmask(:,:,:) * r1_rau0
      DO jk = 2, jpkm1
         znu_w(:,:,jk) = 0.5_wp * ( znu_t(:,:,jk-1) + znu_t(:,:,jk) ) * wmask(:,:,jk)
      END DO

      ! Calculate turbulence intensity parameter Reb
      DO jk = 2, jpkm1
         zReb(:,:,jk) = emix_tmx(:,:,jk) / MAX( 1.e-20_wp, znu_w(:,:,jk) * rn2(:,:,jk) )
      END DO

      ! Define internal wave-induced diffusivity
      DO jk = 2, jpkm1
         zav_wave(:,:,jk) = znu_w(:,:,jk) * zReb(:,:,jk) * r1_6   ! This corresponds to a constant mixing efficiency of 1/6
      END DO

      IF( ln_mevar ) THEN              ! Variable mixing efficiency case : modify zav_wave in the
         DO jk = 2, jpkm1              ! energetic (Reb > 480) and buoyancy-controlled (Reb <10.224 ) regimes
            DO jj = 1, jpj
               DO ji = 1, jpi
                  IF( zReb(ji,jj,jk) > 480.00_wp ) THEN
                     zav_wave(ji,jj,jk) = 3.6515_wp * znu_w(ji,jj,jk) * SQRT( zReb(ji,jj,jk) )
                  ELSEIF( zReb(ji,jj,jk) < 10.224_wp ) THEN
                     zav_wave(ji,jj,jk) = 0.052125_wp * znu_w(ji,jj,jk) * zReb(ji,jj,jk) * SQRT( zReb(ji,jj,jk) )
                  ENDIF
               END DO
            END DO
         END DO
      ENDIF

      DO jk = 2, jpkm1                 ! Bound diffusivity by molecular value and 100 cm2/s
         zav_wave(:,:,jk) = MIN(  MAX( 1.4e-7_wp, zav_wave(:,:,jk) ), 1.e-2_wp  ) * wmask(:,:,jk)
      END DO

      IF( kt == nit000 ) THEN        !* Control print at first time-step: diagnose the energy consumed by zav_wave
         ztpc = 0._wp
         DO jk = 2, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  ztpc = ztpc + fse3w(ji,jj,jk) * e1e2t(ji,jj)   &
                     &         * MAX( 0._wp, rn2(ji,jj,jk) ) * zav_wave(ji,jj,jk) * wmask(ji,jj,jk) * tmask_i(ji,jj)
               END DO
            END DO
         END DO
         IF( lk_mpp )   CALL mpp_sum( ztpc )
         ztpc = rau0 * ztpc ! Global integral of rauo * Kz * N^2 = power contributing to mixing 
 
         IF(lwp) THEN
            WRITE(numout,*)
            WRITE(numout,*) 'zdf_tmx : Internal wave-driven mixing (tmx)'
            WRITE(numout,*) '~~~~~~~ '
            WRITE(numout,*)
            WRITE(numout,*) '      Total power consumption by av_wave: ztpc =  ', ztpc * 1.e-12_wp, 'TW'
         ENDIF
      ENDIF

      !                          ! ----------------------- !
      !                          !   Update  mixing coefs  !                          
      !                          ! ----------------------- !
      !      
      IF( ln_tsdiff ) THEN          !* Option for differential mixing of salinity and temperature
         DO jk = 2, jpkm1              ! Calculate S/T diffusivity ratio as a function of Reb
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zav_ratio(ji,jj,jk) = ( 0.505_wp + 0.495_wp *                                                                  &
                      &   TANH(    0.92_wp * (   LOG10(  MAX( 1.e-20_wp, zReb(ji,jj,jk) * 5._wp * r1_6 )  ) - 0.60_wp   )    )   &
                      &                 ) * wmask(ji,jj,jk)
               END DO
            END DO
         END DO
         CALL iom_put( "av_ratio", zav_ratio )
         DO jk = 2, jpkm1           !* update momentum & tracer diffusivity with wave-driven mixing
            fsavs(:,:,jk) = avt(:,:,jk) + zav_wave(:,:,jk) * zav_ratio(:,:,jk)
            avt  (:,:,jk) = avt(:,:,jk) + zav_wave(:,:,jk)
            avm  (:,:,jk) = avm(:,:,jk) + zav_wave(:,:,jk)
         END DO
         !
      ELSE                          !* update momentum & tracer diffusivity with wave-driven mixing
         DO jk = 2, jpkm1
            fsavs(:,:,jk) = avt(:,:,jk) + zav_wave(:,:,jk)
            avt  (:,:,jk) = avt(:,:,jk) + zav_wave(:,:,jk)
            avm  (:,:,jk) = avm(:,:,jk) + zav_wave(:,:,jk)
         END DO
      ENDIF

      DO jk = 2, jpkm1              !* update momentum diffusivity at wu and wv points
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1  ! vector opt.
               avmu(ji,jj,jk) = avmu(ji,jj,jk) + 0.5_wp * ( zav_wave(ji,jj,jk) + zav_wave(ji+1,jj  ,jk) ) * wumask(ji,jj,jk)
               avmv(ji,jj,jk) = avmv(ji,jj,jk) + 0.5_wp * ( zav_wave(ji,jj,jk) + zav_wave(ji  ,jj+1,jk) ) * wvmask(ji,jj,jk)
            END DO
         END DO
      END DO
      CALL lbc_lnk( avmu, 'U', 1. )   ;   CALL lbc_lnk( avmv, 'V', 1. )      ! lateral boundary condition

      !                             !* output internal wave-driven mixing coefficient
      CALL iom_put( "av_wave", zav_wave )
                                    !* output useful diagnostics: N^2, Kz * N^2 (bflx_tmx), 
                                    !  vertical integral of rau0 * Kz * N^2 (pcmap_tmx), energy density (emix_tmx)
      IF( iom_use("bflx_tmx") .OR. iom_use("pcmap_tmx") ) THEN
         bflx_tmx(:,:,:) = MAX( 0._wp, rn2(:,:,:) ) * zav_wave(:,:,:)
         pcmap_tmx(:,:) = 0._wp
         DO jk = 2, jpkm1
            pcmap_tmx(:,:) = pcmap_tmx(:,:) + fse3w(:,:,jk) * bflx_tmx(:,:,jk) * wmask(:,:,jk)
         END DO
         pcmap_tmx(:,:) = rau0 * pcmap_tmx(:,:)
         CALL iom_put( "bflx_tmx", bflx_tmx )
         CALL iom_put( "pcmap_tmx", pcmap_tmx )
      ENDIF
      CALL iom_put( "bn2", rn2 )
      CALL iom_put( "emix_tmx", emix_tmx )
      
      CALL wrk_dealloc( jpi,jpj,       zfact, zhdep )
      CALL wrk_dealloc( jpi,jpj,jpk,   zwkb, zweight, znu_t, znu_w, zReb )

      IF(ln_ctl)   CALL prt_ctl(tab3d_1=zav_wave , clinfo1=' tmx - av_wave: ', tab3d_2=avt, clinfo2=' avt: ', ovlap=1, kdim=jpk)
      !
      IF( nn_timing == 1 )   CALL timing_stop('zdf_tmx')
      !
   END SUBROUTINE zdf_tmx


   SUBROUTINE zdf_tmx_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE zdf_tmx_init  ***
      !!                     
      !! ** Purpose :   Initialization of the wave-driven vertical mixing, reading
      !!              of input power maps and decay length scales in netcdf files.
      !!
      !! ** Method  : - Read the namzdf_tmx namelist and check the parameters
      !!
      !!              - Read the input data in NetCDF files :
      !!              power available from high-mode wave breaking (mixing_power_bot.nc)
      !!              power available from pycnocline-intensified wave-breaking (mixing_power_pyc.nc)
      !!              power available from critical slope wave-breaking (mixing_power_cri.nc)
      !!              WKB decay scale for high-mode wave-breaking (decay_scale_bot.nc)
      !!              decay scale for critical slope wave-breaking (decay_scale_cri.nc)
      !!
      !! ** input   : - Namlist namzdf_tmx
      !!              - NetCDF files : mixing_power_bot.nc, mixing_power_pyc.nc, mixing_power_cri.nc,
      !!              decay_scale_bot.nc decay_scale_cri.nc
      !!
      !! ** Action  : - Increase by 1 the nstop flag is setting problem encounter
      !!              - Define ebot_tmx, epyc_tmx, ecri_tmx, hbot_tmx, hcri_tmx
      !!
      !! References : de Lavergne et al. 2015, JPO; 2016, in prep.
      !!			  
      !!----------------------------------------------------------------------
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      INTEGER  ::   inum         ! local integer
      INTEGER  ::   ios
      REAL(wp) ::   zbot, zpyc, zcri   ! local scalars
      !!
      NAMELIST/namzdf_tmx_new/ nn_zpyc, ln_mevar, ln_tsdiff
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('zdf_tmx_init')
      !
      REWIND( numnam_ref )              ! Namelist namzdf_tmx in reference namelist : Wave-driven mixing
      READ  ( numnam_ref, namzdf_tmx_new, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namzdf_tmx in reference namelist', lwp )
      !
      REWIND( numnam_cfg )              ! Namelist namzdf_tmx in configuration namelist : Wave-driven mixing
      READ  ( numnam_cfg, namzdf_tmx_new, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namzdf_tmx in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, namzdf_tmx_new )
      !
      IF(lwp) THEN                  ! Control print
         WRITE(numout,*)
         WRITE(numout,*) 'zdf_tmx_init : internal wave-driven mixing'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namzdf_tmx_new : set wave-driven mixing parameters'
         WRITE(numout,*) '      Pycnocline-intensified diss. scales as N (=1) or N^2 (=2) = ', nn_zpyc
         WRITE(numout,*) '      Variable (T) or constant (F) mixing efficiency            = ', ln_mevar
         WRITE(numout,*) '      Differential internal wave-driven mixing (T) or not (F)   = ', ln_tsdiff
      ENDIF
      
      ! The new wave-driven mixing parameterization elevates avt and avm in the interior, and
      ! ensures that avt remains larger than its molecular value (=1.4e-7). Therefore, avtb should 
      ! be set here to a very small value, and avmb to its (uniform) molecular value (=1.4e-6).
      avmb(:) = 1.4e-6_wp        ! viscous molecular value
      avtb(:) = 1.e-10_wp        ! very small diffusive minimum (background avt is specified in zdf_tmx)    
      avtb_2d(:,:) = 1.e0_wp     ! uniform 
      IF(lwp) THEN                  ! Control print
         WRITE(numout,*)
         WRITE(numout,*) '   Force the background value applied to avm & avt in TKE to be everywhere ',   &
            &               'the viscous molecular value & a very small diffusive value, resp.'
      ENDIF
      
      IF( .NOT.lk_zdfddm )   CALL ctl_stop( 'STOP', 'zdf_tmx_init_new : key_zdftmx_new requires key_zdfddm' )
      
      !                             ! allocate tmx arrays
      IF( zdf_tmx_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'zdf_tmx_init : unable to allocate tmx arrays' )
      !
      !                             ! read necessary fields
      CALL iom_open('mixing_power_bot',inum)       ! energy flux for high-mode wave breaking [W/m2]
      CALL iom_get  (inum, jpdom_data, 'field', ebot_tmx, 1 ) 
      CALL iom_close(inum)
      !
      CALL iom_open('mixing_power_pyc',inum)       ! energy flux for pynocline-intensified wave breaking [W/m2]
      CALL iom_get  (inum, jpdom_data, 'field', epyc_tmx, 1 )
      CALL iom_close(inum)
      !
      CALL iom_open('mixing_power_cri',inum)       ! energy flux for critical slope wave breaking [W/m2]
      CALL iom_get  (inum, jpdom_data, 'field', ecri_tmx, 1 )
      CALL iom_close(inum)
      !
      CALL iom_open('decay_scale_bot',inum)        ! spatially variable decay scale for high-mode wave breaking [m]
      CALL iom_get  (inum, jpdom_data, 'field', hbot_tmx, 1 )
      CALL iom_close(inum)
      !
      CALL iom_open('decay_scale_cri',inum)        ! spatially variable decay scale for critical slope wave breaking [m]
      CALL iom_get  (inum, jpdom_data, 'field', hcri_tmx, 1 )
      CALL iom_close(inum)

      ebot_tmx(:,:) = ebot_tmx(:,:) * ssmask(:,:)
      epyc_tmx(:,:) = epyc_tmx(:,:) * ssmask(:,:)
      ecri_tmx(:,:) = ecri_tmx(:,:) * ssmask(:,:)

      ! Set once for all to zero the first and last vertical levels of appropriate variables
      emix_tmx (:,:, 1 ) = 0._wp
      emix_tmx (:,:,jpk) = 0._wp
      zav_ratio(:,:, 1 ) = 0._wp
      zav_ratio(:,:,jpk) = 0._wp
      zav_wave (:,:, 1 ) = 0._wp
      zav_wave (:,:,jpk) = 0._wp

      zbot = glob_sum( e1e2t(:,:) * ebot_tmx(:,:) )
      zpyc = glob_sum( e1e2t(:,:) * epyc_tmx(:,:) )
      zcri = glob_sum( e1e2t(:,:) * ecri_tmx(:,:) )
      IF(lwp) THEN
         WRITE(numout,*) '      High-mode wave-breaking energy:             ', zbot * 1.e-12_wp, 'TW'
         WRITE(numout,*) '      Pycnocline-intensifed wave-breaking energy: ', zpyc * 1.e-12_wp, 'TW'
         WRITE(numout,*) '      Critical slope wave-breaking energy:        ', zcri * 1.e-12_wp, 'TW'
      ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop('zdf_tmx_init')
      !
   END SUBROUTINE zdf_tmx_init

#else
   !!----------------------------------------------------------------------
   !!   Default option          Dummy module                NO Tidal MiXing
   !!----------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_zdftmx = .FALSE.   !: tidal mixing flag
CONTAINS
   SUBROUTINE zdf_tmx_init           ! Dummy routine
      WRITE(*,*) 'zdf_tmx: You should not have seen this print! error?'
   END SUBROUTINE zdf_tmx_init
   SUBROUTINE zdf_tmx( kt )          ! Dummy routine
      WRITE(*,*) 'zdf_tmx: You should not have seen this print! error?', kt
   END SUBROUTINE zdf_tmx
#endif

   !!======================================================================
END MODULE zdftmx
