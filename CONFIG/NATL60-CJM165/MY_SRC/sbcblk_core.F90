MODULE sbcblk_core
   !!======================================================================
   !!                       ***  MODULE  sbcblk_core  ***
   !! Ocean forcing:  momentum, heat and freshwater flux formulation
   !!=====================================================================
   !! History :  1.0  !  2004-08  (U. Schweckendiek)  Original code
   !!            2.0  !  2005-04  (L. Brodeau, A.M. Treguier) additions:
   !!                           -  new bulk routine for efficiency
   !!                           -  WINDS ARE NOW ASSUMED TO BE AT T POINTS in input files !!!!
   !!                           -  file names and file characteristics in namelist
   !!                           -  Implement reading of 6-hourly fields
   !!            3.0  !  2006-06  (G. Madec) sbc rewritting
   !!             -   !  2006-12  (L. Brodeau) Original code for turb_core_2z
   !!            3.2  !  2009-04  (B. Lemaire)  Introduce iom_put
   !!            3.3  !  2010-10  (S. Masson)  add diurnal cycle
   !!            3.4  !  2011-11  (C. Harris) Fill arrays required by CICE
   !!            3.7  !  2014-06  (L. Brodeau) simplification and optimization of CORE bulk
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   sbc_blk_core    : bulk formulation as ocean surface boundary condition (forced mode, CORE bulk formulea)
   !!   blk_oce_core    : computes momentum, heat and freshwater fluxes over ocean
   !!   blk_ice_core    : computes momentum, heat and freshwater fluxes over ice
   !!   turb_core_2z    : Computes turbulent transfert coefficients
   !!   cd_neutral_10m  : Estimate of the neutral drag coefficient at 10m
   !!   psi_m           : universal profile stability function for momentum
   !!   psi_h           : universal profile stability function for temperature and humidity
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE fldread         ! read input fields
   USE sbc_oce         ! Surface boundary condition: ocean fields
   USE cyclone         ! Cyclone 10m wind form trac of cyclone centres
   USE sbcdcy          ! surface boundary condition: diurnal cycle
   USE iom             ! I/O manager library
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! distribued memory computing library
   USE wrk_nemo        ! work arrays
   USE timing          ! Timing
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE prtctl          ! Print control
   USE sbcwave, ONLY   :  cdn_wave ! wave module
   USE sbc_ice         ! Surface boundary condition: ice fields
   USE lib_fortran     ! to use key_nosignedzero
#if defined key_lim3
   USE ice, ONLY       : u_ice, v_ice, jpl, pfrld, a_i_b
   USE limthd_dh       ! for CALL lim_thd_snwblow
#elif defined key_lim2
   USE ice_2, ONLY     : u_ice, v_ice
   USE par_ice_2
#endif

   IMPLICIT NONE
   PRIVATE

   PUBLIC   sbc_blk_core         ! routine called in sbcmod module
#if defined key_lim2 || defined key_lim3
   PUBLIC   blk_ice_core_tau     ! routine called in sbc_ice_lim module
   PUBLIC   blk_ice_core_flx     ! routine called in sbc_ice_lim module
#endif
   PUBLIC   turb_core_2z         ! routine calles in sbcblk_mfs module

! most of the PARAMETERS are made PUBLIC because they are used somehow in the dia_wri routine (flx)      
   INTEGER , PARAMETER         ::   jpfld   = 9   ! maximum number of files to read 
   INTEGER , PUBLIC, PARAMETER ::   jp_wndi = 1   ! index of 10m wind velocity (i-component) (m/s)    at T-point
   INTEGER , PUBLIC, PARAMETER ::   jp_wndj = 2   ! index of 10m wind velocity (j-component) (m/s)    at T-point
   INTEGER , PUBLIC, PARAMETER ::   jp_humi = 3   ! index of specific humidity               ( - )
   INTEGER , PUBLIC, PARAMETER ::   jp_qsr  = 4   ! index of solar heat                      (W/m2)
   INTEGER , PUBLIC, PARAMETER ::   jp_qlw  = 5   ! index of Long wave                       (W/m2)
   INTEGER , PUBLIC, PARAMETER ::   jp_tair = 6   ! index of 10m air temperature             (Kelvin)
   INTEGER , PUBLIC, PARAMETER ::   jp_prec = 7   ! index of total precipitation (rain+snow) (Kg/m2/s)
   INTEGER , PUBLIC, PARAMETER ::   jp_snow = 8   ! index of snow (solid prcipitation)       (kg/m2/s)
   INTEGER , PUBLIC, PARAMETER ::   jp_tdif = 9   ! index of tau diff associated to HF tau   (N/m2)   at T-point

   TYPE(FLD), PUBLIC, ALLOCATABLE, DIMENSION(:) ::   sf   ! structure of input fields (file informations, fields read)
         
   !                                             !!! CORE bulk parameters
   REAL(wp), PARAMETER ::   rhoa =    1.22        ! air density
   REAL(wp), PARAMETER ::   cpa  = 1000.5         ! specific heat of air
   REAL(wp), PARAMETER ::   Lv   =    2.5e6       ! latent heat of vaporization
   REAL(wp), PARAMETER ::   Ls   =    2.839e6     ! latent heat of sublimation
   REAL(wp), PARAMETER ::   Stef =    5.67e-8     ! Stefan Boltzmann constant
   REAL(wp), PARAMETER ::   Cice =    1.4e-3      ! iovi 1.63e-3     ! transfer coefficient over ice
   REAL(wp), PARAMETER ::   albo =    0.066       ! ocean albedo assumed to be constant

   !                                  !!* Namelist namsbc_core : CORE bulk parameters
   LOGICAL  ::   ln_taudif   ! logical flag to use the "mean of stress module - module of mean stress" data
   REAL(wp) ::   rn_pfac     ! multiplication factor for precipitation
   REAL(wp) ::   rn_efac     ! multiplication factor for evaporation (clem)
   REAL(wp) ::   rn_vfac     ! multiplication factor for ice/ocean velocity in the calculation of wind stress (clem)
   REAL(wp) ::   rn_zqt      ! z(q,t) : height of humidity and temperature measurements
   REAL(wp) ::   rn_zu       ! z(u)   : height of wind measurements
   !{ DRAKKAR add JMM
   LOGICAL  ::   ln_kata     ! logical flag for katabatic winds enhancement
   CHARACTER(len=32) :: cl_katfile                !: katabatic filename
   REAL(wp) , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: rmskkatax  !: array for katamask (T-point)
   REAL(wp) , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: rmskkatay  !: array for katamask (T-point)
   !}

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.7 , NEMO-consortium (2014)
   !! $Id: sbcblk_core.F90 5583 2015-07-10 13:23:59Z jchanut $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION sbc_kata_alloc()
      !!----------------------------------------------------------------------
      !!                ***  FUNCTION sbc_kata_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( rmskkatax(jpi,jpj), rmskkatay(jpi,jpj),  STAT= sbc_kata_alloc )
      !
      IF( lk_mpp             )   CALL mpp_sum ( sbc_kata_alloc )
      IF( sbc_kata_alloc > 0 )   CALL ctl_warn('sbc_kata_alloc: allocation of arrays failed')
      !
   END FUNCTION sbc_kata_alloc

   SUBROUTINE sbc_blk_core( kt )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE sbc_blk_core  ***
      !!
      !! ** Purpose :   provide at each time step the surface ocean fluxes
      !!      (momentum, heat, freshwater and runoff)
      !!
      !! ** Method  : (1) READ each fluxes in NetCDF files:
      !!      the 10m wind velocity (i-component) (m/s)    at T-point
      !!      the 10m wind velocity (j-component) (m/s)    at T-point
      !!      the 10m or 2m specific humidity     ( % )
      !!      the solar heat                      (W/m2)
      !!      the Long wave                       (W/m2)
      !!      the 10m or 2m air temperature       (Kelvin)
      !!      the total precipitation (rain+snow) (Kg/m2/s)
      !!      the snow (solid prcipitation)       (kg/m2/s)
      !!      the tau diff associated to HF tau   (N/m2)   at T-point   (ln_taudif=T)
      !!              (2) CALL blk_oce_core
      !!
      !!      C A U T I O N : never mask the surface stress fields
      !!                      the stress is assumed to be in the (i,j) mesh referential
      !!
      !! ** Action  :   defined at each time-step at the air-sea interface
      !!              - utau, vtau  i- and j-component of the wind stress
      !!              - taum        wind stress module at T-point
      !!              - wndm        wind speed  module at T-point over free ocean or leads in presence of sea-ice
      !!              - qns, qsr    non-solar and solar heat fluxes
      !!              - emp         upward mass flux (evapo. - precip.)
      !!              - sfx         salt flux due to freezing/melting (non-zero only if ice is present)
      !!                            (set in limsbc(_2).F90)
      !!
      !! ** References :   Large & Yeager, 2004 / Large & Yeager, 2008
      !!                   Brodeau et al. Ocean Modelling 2010
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time step
      !
      INTEGER  ::   ierror   ! return error code
      INTEGER  ::   ifpr     ! dummy loop indice
      INTEGER  ::   jfld     ! dummy loop arguments
      INTEGER  ::   ios      ! Local integer output status for namelist read
      INTEGER  ::   inum     ! logical unit for IOM
      !!
      CHARACTER(len=100) ::  cn_dir   !   Root directory for location of core files
      TYPE(FLD_N), DIMENSION(jpfld) ::   slf_i     ! array of namelist informations on the fields to read
      TYPE(FLD_N) ::   sn_wndi, sn_wndj, sn_humi, sn_qsr       ! informations about the fields to be read
      TYPE(FLD_N) ::   sn_qlw , sn_tair, sn_prec, sn_snow      !   "                                 "
      TYPE(FLD_N) ::   sn_tdif                                 !   "                                 "
      TYPE(FLD_N) ::   sn_kati, sn_katj                        ! informations for katabatic winds
      NAMELIST/namsbc_core/ cn_dir , ln_taudif, rn_pfac, rn_efac, rn_vfac,  &
         &                  sn_wndi, sn_wndj, sn_humi  , sn_qsr ,           &
         &                  sn_qlw , sn_tair, sn_prec  , sn_snow,           &
         &                  sn_tdif, rn_zqt , rn_zu,                        &
         &                  ln_kata, sn_kati, sn_katj
      !!---------------------------------------------------------------------

      !
      !                                         ! ====================== !
      IF( kt == nit000 ) THEN                   !  First call kt=nit000  !
         !                                      ! ====================== !
         !
         REWIND( numnam_ref )              ! Namelist namsbc_core in reference namelist : CORE bulk parameters
         READ  ( numnam_ref, namsbc_core, IOSTAT = ios, ERR = 901)
901      IF( ios /= 0 ) CALL ctl_nam ( ios , 'namsbc_core in reference namelist', lwp )
         !
         REWIND( numnam_cfg )              ! Namelist namsbc_core in configuration namelist : CORE bulk parameters
         READ  ( numnam_cfg, namsbc_core, IOSTAT = ios, ERR = 902 )
902      IF( ios /= 0 ) CALL ctl_nam ( ios , 'namsbc_core in configuration namelist', lwp )

         IF(lwm) WRITE( numond, namsbc_core )
         !                                         ! check: do we plan to use ln_dm2dc with non-daily forcing?
         IF( ln_dm2dc .AND. sn_qsr%nfreqh /= 24 )   &
            &   CALL ctl_stop( 'sbc_blk_core: ln_dm2dc can be activated only with daily short-wave forcing' )
         IF( ln_dm2dc .AND. sn_qsr%ln_tint ) THEN
            CALL ctl_warn( 'sbc_blk_core: ln_dm2dc is taking care of the temporal interpolation of daily qsr',   &
               &         '              ==> We force time interpolation = .false. for qsr' )
            sn_qsr%ln_tint = .false.
         ENDIF
         !                                         ! store namelist information in an array
         slf_i(jp_wndi) = sn_wndi   ;   slf_i(jp_wndj) = sn_wndj
         slf_i(jp_qsr ) = sn_qsr    ;   slf_i(jp_qlw ) = sn_qlw
         slf_i(jp_tair) = sn_tair   ;   slf_i(jp_humi) = sn_humi
         slf_i(jp_prec) = sn_prec   ;   slf_i(jp_snow) = sn_snow
         slf_i(jp_tdif) = sn_tdif
         !
         lhftau = ln_taudif                        ! do we use HF tau information?
         jfld = jpfld - COUNT( (/.NOT. lhftau/) )
         !
         ALLOCATE( sf(jfld), STAT=ierror )         ! set sf structure
         IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_blk_core: unable to allocate sf structure' )
         DO ifpr= 1, jfld
            ALLOCATE( sf(ifpr)%fnow(jpi,jpj,1) )
            IF( slf_i(ifpr)%ln_tint )   ALLOCATE( sf(ifpr)%fdta(jpi,jpj,1,2) )
         END DO
         !                                         ! fill sf with slf_i and control print
         CALL fld_fill( sf, slf_i, cn_dir, 'sbc_blk_core', 'flux formulation for ocean surface boundary condition', 'namsbc_core' )
         !
         sfx(:,:) = 0._wp                          ! salt flux; zero unless ice is present (computed in limsbc(_2).F90)

         IF ( ln_kata ) THEN
            ! katabatik mask (read in NetCDF file)
            IF ( sbc_kata_alloc() /= 0 ) CALL ctl_stop( 'STOP', 'sbc_kata_alloc: unable to allocate arrays' )
            !
            WRITE(cl_katfile,'(a,a)' ) TRIM( cn_dir ), TRIM( sn_kati%clname )
            CALL iom_open ( cl_katfile, inum )                                       ! open file
            CALL iom_get  ( inum, jpdom_data, sn_kati%clvar, rmskkatax )             ! read the katax array
            CALL iom_get  ( inum, jpdom_data, sn_katj%clvar, rmskkatay )             ! read the katay array
            CALL iom_close( inum )
            ! replace 0 values by 1  in both arrays.
            WHERE (ABS(rmskkatax) < 0.000001 )  ; rmskkatax=1.e0 ; rmskkatay=1.e0 ; END WHERE
            ! check for outliers values
            IF (MAXVAL(rmskkatax) > 6.00001)   CALL ctl_stop( 'Problem  in kata mask : maxval > 6 ')
            IF (MAXVAL(rmskkatay) > 6.00001)   CALL ctl_stop( 'Problem  in kata mask : maxval > 6 ')
         ENDIF
         !
      ENDIF

      CALL fld_read( kt, nn_fsbc, sf )             ! input fields provided at the current time-step

      !                                            ! compute the surface ocean fluxes using CORE bulk formulea
      IF( MOD( kt - 1, nn_fsbc ) == 0 )   CALL blk_oce_core( kt, sf, sst_m, ssu_m, ssv_m )

#if defined key_cice
      IF( MOD( kt - 1, nn_fsbc ) == 0 )   THEN
         qlw_ice(:,:,1)   = sf(jp_qlw)%fnow(:,:,1) 
         qsr_ice(:,:,1)   = sf(jp_qsr)%fnow(:,:,1)
         tatm_ice(:,:)    = sf(jp_tair)%fnow(:,:,1)         
         qatm_ice(:,:)    = sf(jp_humi)%fnow(:,:,1)
         tprecip(:,:)     = sf(jp_prec)%fnow(:,:,1) * rn_pfac
         sprecip(:,:)     = sf(jp_snow)%fnow(:,:,1) * rn_pfac
         wndi_ice(:,:)    = sf(jp_wndi)%fnow(:,:,1)
         wndj_ice(:,:)    = sf(jp_wndj)%fnow(:,:,1)
      ENDIF
#endif
      !
   END SUBROUTINE sbc_blk_core
   
   
   SUBROUTINE blk_oce_core( kt, sf, pst, pu, pv )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE blk_core  ***
      !!
      !! ** Purpose :   provide the momentum, heat and freshwater fluxes at
      !!      the ocean surface at each time step
      !!
      !! ** Method  :   CORE bulk formulea for the ocean using atmospheric
      !!      fields read in sbc_read
      !! 
      !! ** Outputs : - utau    : i-component of the stress at U-point  (N/m2)
      !!              - vtau    : j-component of the stress at V-point  (N/m2)
      !!              - taum    : Wind stress module at T-point         (N/m2)
      !!              - wndm    : Wind speed module at T-point          (m/s)
      !!              - qsr     : Solar heat flux over the ocean        (W/m2)
      !!              - qns     : Non Solar heat flux over the ocean    (W/m2)
      !!              - emp     : evaporation minus precipitation       (kg/m2/s)
      !!
      !!  ** Nota  :   sf has to be a dummy argument for AGRIF on NEC
      !!---------------------------------------------------------------------
      INTEGER  , INTENT(in   )                 ::   kt    ! time step index
      TYPE(fld), INTENT(inout), DIMENSION(:)   ::   sf    ! input data
      REAL(wp) , INTENT(in)   , DIMENSION(:,:) ::   pst   ! surface temperature                      [Celcius]
      REAL(wp) , INTENT(in)   , DIMENSION(:,:) ::   pu    ! surface current at U-point (i-component) [m/s]
      REAL(wp) , INTENT(in)   , DIMENSION(:,:) ::   pv    ! surface current at V-point (j-component) [m/s]
      !
      INTEGER  ::   ji, jj               ! dummy loop indices
      REAL(wp) ::   zcoef_qsatw, zztmp   ! local variable
      REAL(wp), DIMENSION(:,:), POINTER ::   zwnd_i, zwnd_j    ! wind speed components at T-point
      REAL(wp), DIMENSION(:,:), POINTER ::   zqsatw            ! specific humidity at pst
      REAL(wp), DIMENSION(:,:), POINTER ::   zqlw, zqsb        ! long wave and sensible heat fluxes
      REAL(wp), DIMENSION(:,:), POINTER ::   zqla, zevap       ! latent heat fluxes and evaporation
      REAL(wp), DIMENSION(:,:), POINTER ::   Cd                ! transfer coefficient for momentum      (tau)
      REAL(wp), DIMENSION(:,:), POINTER ::   Ch                ! transfer coefficient for sensible heat (Q_sens)
      REAL(wp), DIMENSION(:,:), POINTER ::   Ce                ! tansfert coefficient for evaporation   (Q_lat)
      REAL(wp), DIMENSION(:,:), POINTER ::   zst               ! surface temperature in Kelvin
      REAL(wp), DIMENSION(:,:), POINTER ::   zt_zu             ! air temperature at wind speed height
      REAL(wp), DIMENSION(:,:), POINTER ::   zq_zu             ! air spec. hum.  at wind speed height
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('blk_oce_core')
      !
      CALL wrk_alloc( jpi,jpj, zwnd_i, zwnd_j, zqsatw, zqlw, zqsb, zqla, zevap )
      CALL wrk_alloc( jpi,jpj, Cd, Ch, Ce, zst, zt_zu, zq_zu )
      !
      ! local scalars ( place there for vector optimisation purposes)
      zcoef_qsatw = 0.98 * 640380. / rhoa
      
      zst(:,:) = pst(:,:) + rt0      ! convert SST from  Celsius to Kelvin (and set minimum value far above 0 K)

      ! ----------------------------------------------------------------------------- !
      !      0   Wind components and module at T-point relative to the moving ocean   !
      ! ----------------------------------------------------------------------------- !

      ! ... components ( U10m - U_oce ) at T-point (unmasked)
      zwnd_i(:,:) = 0.e0  
      zwnd_j(:,:) = 0.e0
#if defined key_cyclone
      CALL wnd_cyc( kt, zwnd_i, zwnd_j )    ! add analytical tropical cyclone (Vincent et al. JGR 2012)
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! vect. opt.
            sf(jp_wndi)%fnow(ji,jj,1) = sf(jp_wndi)%fnow(ji,jj,1) + zwnd_i(ji,jj)
            sf(jp_wndj)%fnow(ji,jj,1) = sf(jp_wndj)%fnow(ji,jj,1) + zwnd_j(ji,jj)
         END DO
      END DO
#endif
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! vect. opt.
            zwnd_i(ji,jj) = (  sf(jp_wndi)%fnow(ji,jj,1) - rn_vfac * 0.5 * ( pu(ji-1,jj  ) + pu(ji,jj) )  )
            zwnd_j(ji,jj) = (  sf(jp_wndj)%fnow(ji,jj,1) - rn_vfac * 0.5 * ( pv(ji  ,jj-1) + pv(ji,jj) )  )
         END DO
      END DO
      CALL lbc_lnk( zwnd_i(:,:) , 'T', -1. )
      CALL lbc_lnk( zwnd_j(:,:) , 'T', -1. )
      ! ... scalar wind ( = | U10m - U_oce | ) at T-point (masked)
      wndm(:,:) = SQRT(  zwnd_i(:,:) * zwnd_i(:,:)   &
         &             + zwnd_j(:,:) * zwnd_j(:,:)  ) * tmask(:,:,1)

      ! ----------------------------------------------------------------------------- !
      !      I   Radiative FLUXES                                                     !
      ! ----------------------------------------------------------------------------- !

      ! ocean albedo assumed to be constant + modify now Qsr to include the diurnal cycle                    ! Short Wave
      zztmp = 1. - albo
      IF( ln_dm2dc ) THEN   ;   qsr(:,:) = zztmp * sbc_dcy( sf(jp_qsr)%fnow(:,:,1) ) * tmask(:,:,1)
      ELSE                  ;   qsr(:,:) = zztmp *          sf(jp_qsr)%fnow(:,:,1)   * tmask(:,:,1)
      ENDIF

      zqlw(:,:) = (  sf(jp_qlw)%fnow(:,:,1) - Stef * zst(:,:)*zst(:,:)*zst(:,:)*zst(:,:)  ) * tmask(:,:,1)   ! Long  Wave
      ! ----------------------------------------------------------------------------- !
      !     II    Turbulent FLUXES                                                    !
      ! ----------------------------------------------------------------------------- !

      ! ... specific humidity at SST and IST
      zqsatw(:,:) = zcoef_qsatw * EXP( -5107.4 / zst(:,:) )

      ! ... NCAR Bulk formulae, computation of Cd, Ch, Ce at T-point :
      CALL turb_core_2z( rn_zqt, rn_zu, zst, sf(jp_tair)%fnow, zqsatw, sf(jp_humi)%fnow, wndm,   &
         &               Cd, Ch, Ce, zt_zu, zq_zu )
    
      ! ... tau module, i and j component
      DO jj = 1, jpj
         DO ji = 1, jpi
            zztmp = rhoa * wndm(ji,jj) * Cd(ji,jj)
            taum  (ji,jj) = zztmp * wndm  (ji,jj)
            zwnd_i(ji,jj) = zztmp * zwnd_i(ji,jj)
            zwnd_j(ji,jj) = zztmp * zwnd_j(ji,jj)
         END DO
      END DO

      IF ( ln_kata )   THEN   ! correction of wind stress at T point only ( no impact on wind speed)
        zwnd_i(:,:) = rmskkatax(:,:)*zwnd_i(:,:)
        zwnd_j(:,:) = rmskkatay(:,:)*zwnd_j(:,:)
      ENDIF
      ! ... add the HF tau contribution to the wind stress module?
      IF( lhftau ) THEN
         taum(:,:) = taum(:,:) + sf(jp_tdif)%fnow(:,:,1)
      ENDIF
      CALL iom_put( "taum_oce", taum )   ! output wind stress module

      ! ... utau, vtau at U- and V_points, resp.
      !     Note the use of 0.5*(2-umask) in order to unmask the stress along coastlines
      !     Note the use of MAX(tmask(i,j),tmask(i+1,j) is to mask tau over ice shelves
      DO jj = 1, jpjm1
         DO ji = 1, fs_jpim1
            utau(ji,jj) = 0.5 * ( 2. - umask(ji,jj,1) ) * ( zwnd_i(ji,jj) + zwnd_i(ji+1,jj  ) ) &
               &          * MAX(tmask(ji,jj,1),tmask(ji+1,jj,1))
            vtau(ji,jj) = 0.5 * ( 2. - vmask(ji,jj,1) ) * ( zwnd_j(ji,jj) + zwnd_j(ji  ,jj+1) ) &
               &          * MAX(tmask(ji,jj,1),tmask(ji,jj+1,1))
         END DO
      END DO
      CALL lbc_lnk( utau(:,:), 'U', -1. )
      CALL lbc_lnk( vtau(:,:), 'V', -1. )

    
      !  Turbulent fluxes over ocean
      ! -----------------------------
      IF( ABS( rn_zu - rn_zqt) < 0.01_wp ) THEN
         !! q_air and t_air are (or "are almost") given at 10m (wind reference height)
         zevap(:,:) = rn_efac*MAX( 0._wp,     rhoa*Ce(:,:)*( zqsatw(:,:) - sf(jp_humi)%fnow(:,:,1) )*wndm(:,:) ) ! Evaporation
         zqsb (:,:) =                     cpa*rhoa*Ch(:,:)*( zst   (:,:) - sf(jp_tair)%fnow(:,:,1) )*wndm(:,:)   ! Sensible Heat
      ELSE
         !! q_air and t_air are not given at 10m (wind reference height)
         ! Values of temp. and hum. adjusted to height of wind during bulk algorithm iteration must be used!!!
         zevap(:,:) = rn_efac*MAX( 0._wp,     rhoa*Ce(:,:)*( zqsatw(:,:) - zq_zu(:,:) )*wndm(:,:) )   ! Evaporation
         zqsb (:,:) =                     cpa*rhoa*Ch(:,:)*( zst   (:,:) - zt_zu(:,:) )*wndm(:,:)     ! Sensible Heat
      ENDIF
      zqla (:,:) = Lv * zevap(:,:)                                                              ! Latent Heat
      IF(ln_ctl) THEN
         CALL prt_ctl( tab2d_1=zqla  , clinfo1=' blk_oce_core: zqla   : ', tab2d_2=Ce , clinfo2=' Ce  : ' )
         CALL prt_ctl( tab2d_1=zqsb  , clinfo1=' blk_oce_core: zqsb   : ', tab2d_2=Ch , clinfo2=' Ch  : ' )
         CALL prt_ctl( tab2d_1=zqlw  , clinfo1=' blk_oce_core: zqlw   : ', tab2d_2=qsr, clinfo2=' qsr : ' )
         CALL prt_ctl( tab2d_1=zqsatw, clinfo1=' blk_oce_core: zqsatw : ', tab2d_2=zst, clinfo2=' zst : ' )
         CALL prt_ctl( tab2d_1=utau  , clinfo1=' blk_oce_core: utau   : ', mask1=umask,   &
            &          tab2d_2=vtau  , clinfo2=              ' vtau : '  , mask2=vmask )
         CALL prt_ctl( tab2d_1=wndm  , clinfo1=' blk_oce_core: wndm   : ')
         CALL prt_ctl( tab2d_1=zst   , clinfo1=' blk_oce_core: zst    : ')
      ENDIF
       
      ! ----------------------------------------------------------------------------- !
      !     III    Total FLUXES                                                       !
      ! ----------------------------------------------------------------------------- !
      !
      emp (:,:) = (  zevap(:,:)                                          &   ! mass flux (evap. - precip.)
         &         - sf(jp_prec)%fnow(:,:,1) * rn_pfac  ) * tmask(:,:,1)
      !
      qns(:,:) = zqlw(:,:) - zqsb(:,:) - zqla(:,:)                                &   ! Downward Non Solar 
         &     - sf(jp_snow)%fnow(:,:,1) * rn_pfac * lfus                         &   ! remove latent melting heat for solid precip
         &     - zevap(:,:) * pst(:,:) * rcp                                      &   ! remove evap heat content at SST
         &     + ( sf(jp_prec)%fnow(:,:,1) - sf(jp_snow)%fnow(:,:,1) ) * rn_pfac  &   ! add liquid precip heat content at Tair
         &     * ( sf(jp_tair)%fnow(:,:,1) - rt0 ) * rcp                          &
         &     + sf(jp_snow)%fnow(:,:,1) * rn_pfac                                &   ! add solid  precip heat content at min(Tair,Tsnow)
         &     * ( MIN( sf(jp_tair)%fnow(:,:,1), rt0_snow ) - rt0 ) * cpic * tmask(:,:,1)
      !
#if defined key_lim3
      qns_oce(:,:) = zqlw(:,:) - zqsb(:,:) - zqla(:,:)                                ! non solar without emp (only needed by LIM3)
      qsr_oce(:,:) = qsr(:,:)
#endif
      !
      IF ( nn_ice == 0 ) THEN
         CALL iom_put( "qlw_oce" ,   zqlw )                 ! output downward longwave heat over the ocean
         CALL iom_put( "qsb_oce" , - zqsb )                 ! output downward sensible heat over the ocean
         CALL iom_put( "qla_oce" , - zqla )                 ! output downward latent   heat over the ocean
         CALL iom_put( "qemp_oce",   qns-zqlw+zqsb+zqla )   ! output downward heat content of E-P over the ocean
         CALL iom_put( "qns_oce" ,   qns  )                 ! output downward non solar heat over the ocean
         CALL iom_put( "qsr_oce" ,   qsr  )                 ! output downward solar heat over the ocean
         CALL iom_put( "qt_oce"  ,   qns+qsr )              ! output total downward heat over the ocean
         tprecip(:,:) = sf(jp_prec)%fnow(:,:,1) * rn_pfac   ! output total precipitation [kg/m2/s]
         sprecip(:,:) = sf(jp_snow)%fnow(:,:,1) * rn_pfac   ! output solid precipitation [kg/m2/s]
         CALL iom_put( 'snowpre', sprecip * 86400. )        ! Snow
         CALL iom_put( 'precip' , tprecip * 86400. )        ! Total precipitation
      ENDIF
      !
      IF(ln_ctl) THEN
         CALL prt_ctl(tab2d_1=zqsb , clinfo1=' blk_oce_core: zqsb   : ', tab2d_2=zqlw , clinfo2=' zqlw  : ')
         CALL prt_ctl(tab2d_1=zqla , clinfo1=' blk_oce_core: zqla   : ', tab2d_2=qsr  , clinfo2=' qsr   : ')
         CALL prt_ctl(tab2d_1=pst  , clinfo1=' blk_oce_core: pst    : ', tab2d_2=emp  , clinfo2=' emp   : ')
         CALL prt_ctl(tab2d_1=utau , clinfo1=' blk_oce_core: utau   : ', mask1=umask,   &
            &         tab2d_2=vtau , clinfo2=              ' vtau  : ' , mask2=vmask )
      ENDIF
      !
      CALL wrk_dealloc( jpi,jpj, zwnd_i, zwnd_j, zqsatw, zqlw, zqsb, zqla, zevap )
      CALL wrk_dealloc( jpi,jpj, Cd, Ch, Ce, zst, zt_zu, zq_zu )
      !
      IF( nn_timing == 1 )  CALL timing_stop('blk_oce_core')
      !
   END SUBROUTINE blk_oce_core
 
   
#if defined key_lim2 || defined key_lim3
   SUBROUTINE blk_ice_core_tau
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE blk_ice_core_tau  ***
      !!
      !! ** Purpose :   provide the surface boundary condition over sea-ice
      !!
      !! ** Method  :   compute momentum using CORE bulk
      !!                formulea, ice variables and read atmospheric fields.
      !!                NB: ice drag coefficient is assumed to be a constant
      !!---------------------------------------------------------------------
      INTEGER  ::   ji, jj    ! dummy loop indices
      REAL(wp) ::   zcoef_wnorm, zcoef_wnorm2
      REAL(wp) ::   zwnorm_f, zwndi_f , zwndj_f               ! relative wind module and components at F-point
      REAL(wp) ::             zwndi_t , zwndj_t               ! relative wind components at T-point
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('blk_ice_core_tau')
      !
      ! local scalars ( place there for vector optimisation purposes)
      zcoef_wnorm  = rhoa * Cice
      zcoef_wnorm2 = rhoa * Cice * 0.5

!!gm brutal....
      utau_ice  (:,:) = 0._wp
      vtau_ice  (:,:) = 0._wp
      wndm_ice  (:,:) = 0._wp
!!gm end

      ! ----------------------------------------------------------------------------- !
      !    Wind components and module relative to the moving ocean ( U10m - U_ice )   !
      ! ----------------------------------------------------------------------------- !
      SELECT CASE( cp_ice_msh )
      CASE( 'I' )                  ! B-grid ice dynamics :   I-point (i.e. F-point with sea-ice indexation)
         !                           and scalar wind at T-point ( = | U10m - U_ice | ) (masked)
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! B grid : NO vector opt
               ! ... scalar wind at I-point (fld being at T-point)
               zwndi_f = 0.25 * (  sf(jp_wndi)%fnow(ji-1,jj  ,1) + sf(jp_wndi)%fnow(ji  ,jj  ,1)   &
                  &              + sf(jp_wndi)%fnow(ji-1,jj-1,1) + sf(jp_wndi)%fnow(ji  ,jj-1,1)  ) - rn_vfac * u_ice(ji,jj)
               zwndj_f = 0.25 * (  sf(jp_wndj)%fnow(ji-1,jj  ,1) + sf(jp_wndj)%fnow(ji  ,jj  ,1)   &
                  &              + sf(jp_wndj)%fnow(ji-1,jj-1,1) + sf(jp_wndj)%fnow(ji  ,jj-1,1)  ) - rn_vfac * v_ice(ji,jj)
               zwnorm_f = zcoef_wnorm * SQRT( zwndi_f * zwndi_f + zwndj_f * zwndj_f )
               ! ... ice stress at I-point
               utau_ice(ji,jj) = zwnorm_f * zwndi_f
               vtau_ice(ji,jj) = zwnorm_f * zwndj_f
               IF (ln_kata ) THEN   ! Apply Katabatic wind enhancement on stress
                  utau_ice(ji,jj) = utau_ice(ji,jj) * 0.25 * ( rmskkatax(ji-1,jj  ) + rmskkatax(ji,jj  )  &
                     &                                   + rmskkatax(ji-1,jj-1) + rmskkatax(ji,jj-1) )
                  vtau_ice(ji,jj) = vtau_ice(ji,jj) * 0.25 * ( rmskkatay(ji-1,jj  ) + rmskkatay(ji,jj  )  &
                     &                                   + rmskkatay(ji-1,jj-1) + rmskkatay(ji,jj-1) )
               ENDIF
               ! ... scalar wind at T-point (fld being at T-point)
               zwndi_t = sf(jp_wndi)%fnow(ji,jj,1) - rn_vfac * 0.25 * (  u_ice(ji,jj+1) + u_ice(ji+1,jj+1)   &
                  &                                                    + u_ice(ji,jj  ) + u_ice(ji+1,jj  )  )
               zwndj_t = sf(jp_wndj)%fnow(ji,jj,1) - rn_vfac * 0.25 * (  v_ice(ji,jj+1) + v_ice(ji+1,jj+1)   &
                  &                                                    + v_ice(ji,jj  ) + v_ice(ji+1,jj  )  )
               wndm_ice(ji,jj)  = SQRT( zwndi_t * zwndi_t + zwndj_t * zwndj_t ) * tmask(ji,jj,1)
            END DO
         END DO
         CALL lbc_lnk( utau_ice, 'I', -1. )
         CALL lbc_lnk( vtau_ice, 'I', -1. )
         CALL lbc_lnk( wndm_ice, 'T',  1. )
         !
      CASE( 'C' )                  ! C-grid ice dynamics :   U & V-points (same as ocean)
         DO jj = 2, jpj
            DO ji = fs_2, jpi   ! vect. opt.
               zwndi_t = (  sf(jp_wndi)%fnow(ji,jj,1) - rn_vfac * 0.5 * ( u_ice(ji-1,jj  ) + u_ice(ji,jj) )  )
               zwndj_t = (  sf(jp_wndj)%fnow(ji,jj,1) - rn_vfac * 0.5 * ( v_ice(ji  ,jj-1) + v_ice(ji,jj) )  )
               wndm_ice(ji,jj) = SQRT( zwndi_t * zwndi_t + zwndj_t * zwndj_t ) * tmask(ji,jj,1)
            END DO
         END DO
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vect. opt.
               utau_ice(ji,jj) = zcoef_wnorm2 * ( wndm_ice(ji+1,jj  ) + wndm_ice(ji,jj) )                          &
                  &          * ( 0.5 * (sf(jp_wndi)%fnow(ji+1,jj,1) + sf(jp_wndi)%fnow(ji,jj,1) ) - rn_vfac * u_ice(ji,jj) )
               vtau_ice(ji,jj) = zcoef_wnorm2 * ( wndm_ice(ji,jj+1  ) + wndm_ice(ji,jj) )                          &
                  &          * ( 0.5 * (sf(jp_wndj)%fnow(ji,jj+1,1) + sf(jp_wndj)%fnow(ji,jj,1) ) - rn_vfac * v_ice(ji,jj) )
               IF (ln_kata ) THEN ! Apply Katabatic wind enhancement on stress
                  utau_ice(ji,jj) = utau_ice(ji,jj) * 0.5 * ( rmskkatax(ji+1,jj) + rmskkatax(ji,jj) )
                  vtau_ice(ji,jj) = vtau_ice(ji,jj) * 0.5 * ( rmskkatay(ji,jj+1) + rmskkatay(ji,jj) )
               ENDIF
            END DO
         END DO
         CALL lbc_lnk( utau_ice, 'U', -1. )
         CALL lbc_lnk( vtau_ice, 'V', -1. )
         CALL lbc_lnk( wndm_ice, 'T',  1. )
         !
      END SELECT

      IF(ln_ctl) THEN
         CALL prt_ctl(tab2d_1=utau_ice  , clinfo1=' blk_ice_core: utau_ice : ', tab2d_2=vtau_ice  , clinfo2=' vtau_ice : ')
         CALL prt_ctl(tab2d_1=wndm_ice  , clinfo1=' blk_ice_core: wndm_ice : ')
      ENDIF

      IF( nn_timing == 1 )  CALL timing_stop('blk_ice_core_tau')
      
   END SUBROUTINE blk_ice_core_tau


   SUBROUTINE blk_ice_core_flx( ptsu, palb )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE blk_ice_core_flx  ***
      !!
      !! ** Purpose :   provide the surface boundary condition over sea-ice
      !!
      !! ** Method  :   compute heat and freshwater exchanged
      !!                between atmosphere and sea-ice using CORE bulk
      !!                formulea, ice variables and read atmmospheric fields.
      !! 
      !! caution : the net upward water flux has with mm/day unit
      !!---------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:,:), INTENT(in)  ::   ptsu          ! sea ice surface temperature
      REAL(wp), DIMENSION(:,:,:), INTENT(in)  ::   palb          ! ice albedo (all skies)
      !!
      INTEGER  ::   ji, jj, jl    ! dummy loop indices
      REAL(wp) ::   zst2, zst3
      REAL(wp) ::   zcoef_dqlw, zcoef_dqla, zcoef_dqsb
      REAL(wp) ::   zztmp, z1_lsub                               ! temporary variable
      !!
      REAL(wp), DIMENSION(:,:,:), POINTER ::   z_qlw             ! long wave heat flux over ice
      REAL(wp), DIMENSION(:,:,:), POINTER ::   z_qsb             ! sensible  heat flux over ice
      REAL(wp), DIMENSION(:,:,:), POINTER ::   z_dqlw            ! long wave heat sensitivity over ice
      REAL(wp), DIMENSION(:,:,:), POINTER ::   z_dqsb            ! sensible  heat sensitivity over ice
      REAL(wp), DIMENSION(:,:)  , POINTER ::   zevap, zsnw       ! evaporation and snw distribution after wind blowing (LIM3)
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('blk_ice_core_flx')
      !
      CALL wrk_alloc( jpi,jpj,jpl, z_qlw, z_qsb, z_dqlw, z_dqsb ) 

      ! local scalars ( place there for vector optimisation purposes)
      zcoef_dqlw   = 4.0 * 0.95 * Stef
      zcoef_dqla   = -Ls * Cice * 11637800. * (-5897.8)
      zcoef_dqsb   = rhoa * cpa * Cice

      zztmp = 1. / ( 1. - albo )
      !                                     ! ========================== !
      DO jl = 1, jpl                        !  Loop over ice categories  !
         !                                  ! ========================== !
         DO jj = 1 , jpj
            DO ji = 1, jpi
               ! ----------------------------!
               !      I   Radiative FLUXES   !
               ! ----------------------------!
               zst2 = ptsu(ji,jj,jl) * ptsu(ji,jj,jl)
               zst3 = ptsu(ji,jj,jl) * zst2
               ! Short Wave (sw)
               qsr_ice(ji,jj,jl) = zztmp * ( 1. - palb(ji,jj,jl) ) * qsr(ji,jj)
               ! Long  Wave (lw)
               !  0.95 is the emissivity 
               z_qlw(ji,jj,jl) = 0.95 * ( sf(jp_qlw)%fnow(ji,jj,1) - Stef * ptsu(ji,jj,jl) * zst3 ) * tmask(ji,jj,1)
               ! lw sensitivity
               z_dqlw(ji,jj,jl) = zcoef_dqlw * zst3                                               

               ! ----------------------------!
               !     II    Turbulent FLUXES  !
               ! ----------------------------!

               ! ... turbulent heat fluxes
               ! Sensible Heat
               z_qsb(ji,jj,jl) = rhoa * cpa * Cice * wndm_ice(ji,jj) * ( ptsu(ji,jj,jl) - sf(jp_tair)%fnow(ji,jj,1) )
               ! Latent Heat
               qla_ice(ji,jj,jl) = rn_efac * MAX( 0.e0, rhoa * Ls  * Cice * wndm_ice(ji,jj)   &                           
                  &                         * (  11637800. * EXP( -5897.8 / ptsu(ji,jj,jl) ) / rhoa - sf(jp_humi)%fnow(ji,jj,1)  ) )
              ! Latent heat sensitivity for ice (Dqla/Dt)
               IF( qla_ice(ji,jj,jl) > 0._wp ) THEN
                  dqla_ice(ji,jj,jl) = rn_efac * zcoef_dqla * wndm_ice(ji,jj) / ( zst2 ) * EXP( -5897.8 / ptsu(ji,jj,jl) )
               ELSE
                  dqla_ice(ji,jj,jl) = 0._wp
               ENDIF

               ! Sensible heat sensitivity (Dqsb_ice/Dtn_ice)
               z_dqsb(ji,jj,jl) = zcoef_dqsb * wndm_ice(ji,jj)

               ! ----------------------------!
               !     III    Total FLUXES     !
               ! ----------------------------!
               ! Downward Non Solar flux
               qns_ice (ji,jj,jl) =     z_qlw (ji,jj,jl) - z_qsb (ji,jj,jl) - qla_ice (ji,jj,jl)
               ! Total non solar heat flux sensitivity for ice
               dqns_ice(ji,jj,jl) = - ( z_dqlw(ji,jj,jl) + z_dqsb(ji,jj,jl) + dqla_ice(ji,jj,jl) )
            END DO
            !
         END DO
         !
      END DO
      !
      tprecip(:,:) = sf(jp_prec)%fnow(:,:,1) * rn_pfac      ! total precipitation [kg/m2/s]
      sprecip(:,:) = sf(jp_snow)%fnow(:,:,1) * rn_pfac      ! solid precipitation [kg/m2/s]
      CALL iom_put( 'snowpre', sprecip * 86400. )                  ! Snow precipitation
      CALL iom_put( 'precip' , tprecip * 86400. )                  ! Total precipitation

#if defined  key_lim3
      CALL wrk_alloc( jpi,jpj, zevap, zsnw ) 

      ! --- evaporation --- !
      z1_lsub = 1._wp / Lsub
      evap_ice (:,:,:) = qla_ice (:,:,:) * z1_lsub ! sublimation
      devap_ice(:,:,:) = dqla_ice(:,:,:) * z1_lsub
      zevap    (:,:)   = emp(:,:) + tprecip(:,:)   ! evaporation over ocean

      ! --- evaporation minus precipitation --- !
      zsnw(:,:) = 0._wp
      CALL lim_thd_snwblow( pfrld, zsnw )  ! snow distribution over ice after wind blowing 
      emp_oce(:,:) = pfrld(:,:) * zevap(:,:) - ( tprecip(:,:) - sprecip(:,:) ) - sprecip(:,:) * (1._wp - zsnw )
      emp_ice(:,:) = SUM( a_i_b(:,:,:) * evap_ice(:,:,:), dim=3 ) - sprecip(:,:) * zsnw
      emp_tot(:,:) = emp_oce(:,:) + emp_ice(:,:)

      ! --- heat flux associated with emp --- !
      qemp_oce(:,:) = - pfrld(:,:) * zevap(:,:) * sst_m(:,:) * rcp                               & ! evap at sst
         &          + ( tprecip(:,:) - sprecip(:,:) ) * ( sf(jp_tair)%fnow(:,:,1) - rt0 ) * rcp  & ! liquid precip at Tair
         &          +   sprecip(:,:) * ( 1._wp - zsnw ) *                                        & ! solid precip at min(Tair,Tsnow)
         &              ( ( MIN( sf(jp_tair)%fnow(:,:,1), rt0_snow ) - rt0 ) * cpic * tmask(:,:,1) - lfus )
      qemp_ice(:,:) =   sprecip(:,:) * zsnw *                                                    & ! solid precip (only)
         &              ( ( MIN( sf(jp_tair)%fnow(:,:,1), rt0_snow ) - rt0 ) * cpic * tmask(:,:,1) - lfus )

      ! --- total solar and non solar fluxes --- !
      qns_tot(:,:) = pfrld(:,:) * qns_oce(:,:) + SUM( a_i_b(:,:,:) * qns_ice(:,:,:), dim=3 ) + qemp_ice(:,:) + qemp_oce(:,:)
      qsr_tot(:,:) = pfrld(:,:) * qsr_oce(:,:) + SUM( a_i_b(:,:,:) * qsr_ice(:,:,:), dim=3 )

      ! --- heat content of precip over ice in J/m3 (to be used in 1D-thermo) --- !
      qprec_ice(:,:) = rhosn * ( ( MIN( sf(jp_tair)%fnow(:,:,1), rt0_snow ) - rt0 ) * cpic * tmask(:,:,1) - lfus )

      CALL wrk_dealloc( jpi,jpj, zevap, zsnw ) 
#endif

      !--------------------------------------------------------------------
      ! FRACTIONs of net shortwave radiation which is not absorbed in the
      ! thin surface layer and penetrates inside the ice cover
      ! ( Maykut and Untersteiner, 1971 ; Ebert and Curry, 1993 )
      !
      fr1_i0(:,:) = ( 0.18 * ( 1.0 - cldf_ice ) + 0.35 * cldf_ice )
      fr2_i0(:,:) = ( 0.82 * ( 1.0 - cldf_ice ) + 0.65 * cldf_ice )
      !
      !
      IF(ln_ctl) THEN
         CALL prt_ctl(tab3d_1=qla_ice , clinfo1=' blk_ice_core: qla_ice  : ', tab3d_2=z_qsb   , clinfo2=' z_qsb    : ', kdim=jpl)
         CALL prt_ctl(tab3d_1=z_qlw   , clinfo1=' blk_ice_core: z_qlw    : ', tab3d_2=dqla_ice, clinfo2=' dqla_ice : ', kdim=jpl)
         CALL prt_ctl(tab3d_1=z_dqsb  , clinfo1=' blk_ice_core: z_dqsb   : ', tab3d_2=z_dqlw  , clinfo2=' z_dqlw   : ', kdim=jpl)
         CALL prt_ctl(tab3d_1=dqns_ice, clinfo1=' blk_ice_core: dqns_ice : ', tab3d_2=qsr_ice , clinfo2=' qsr_ice  : ', kdim=jpl)
         CALL prt_ctl(tab3d_1=ptsu    , clinfo1=' blk_ice_core: ptsu     : ', tab3d_2=qns_ice , clinfo2=' qns_ice  : ', kdim=jpl)
         CALL prt_ctl(tab2d_1=tprecip , clinfo1=' blk_ice_core: tprecip  : ', tab2d_2=sprecip , clinfo2=' sprecip  : ')
      ENDIF

      CALL wrk_dealloc( jpi,jpj,jpl, z_qlw, z_qsb, z_dqlw, z_dqsb )
      !
      IF( nn_timing == 1 )  CALL timing_stop('blk_ice_core_flx')
      
   END SUBROUTINE blk_ice_core_flx
#endif

   SUBROUTINE turb_core_2z( zt, zu, sst, T_zt, q_sat, q_zt, dU,    &
      &                      Cd, Ch, Ce , T_zu, q_zu )
      !!----------------------------------------------------------------------
      !!                      ***  ROUTINE  turb_core  ***
      !!
      !! ** Purpose :   Computes turbulent transfert coefficients of surface
      !!                fluxes according to Large & Yeager (2004) and Large & Yeager (2008)
      !!                If relevant (zt /= zu), adjust temperature and humidity from height zt to zu
      !!
      !! ** Method : Monin Obukhov Similarity Theory 
      !!             + Large & Yeager (2004,2008) closure: CD_n10 = f(U_n10)
      !!
      !! ** References :   Large & Yeager, 2004 / Large & Yeager, 2008
      !!
      !! ** Last update: Laurent Brodeau, June 2014:
      !!    - handles both cases zt=zu and zt/=zu
      !!    - optimized: less 2D arrays allocated and less operations
      !!    - better first guess of stability by checking air-sea difference of virtual temperature
      !!       rather than temperature difference only...
      !!    - added function "cd_neutral_10m" that uses the improved parametrization of 
      !!      Large & Yeager 2008. Drag-coefficient reduction for Cyclone conditions!
      !!    - using code-wide physical constants defined into "phycst.mod" rather than redifining them
      !!      => 'vkarmn' and 'grav'
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(in   )                     ::   zt       ! height for T_zt and q_zt                   [m]
      REAL(wp), INTENT(in   )                     ::   zu       ! height for dU                              [m]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   sst      ! sea surface temperature              [Kelvin]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   T_zt     ! potential air temperature            [Kelvin]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   q_sat    ! sea surface specific humidity         [kg/kg]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   q_zt     ! specific air humidity                 [kg/kg]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   dU       ! relative wind module at zu            [m/s]
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   Cd       ! transfer coefficient for momentum         (tau)
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   Ch       ! transfer coefficient for sensible heat (Q_sens)
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   Ce       ! transfert coefficient for evaporation   (Q_lat)
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   T_zu     ! air temp. shifted at zu                     [K]
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   q_zu     ! spec. hum.  shifted at zu               [kg/kg]
      !
      INTEGER ::   j_itt
      INTEGER , PARAMETER ::   nb_itt = 5       ! number of itterations
      LOGICAL ::   l_zt_equal_zu = .FALSE.      ! if q and t are given at different height than U
      !
      REAL(wp), DIMENSION(:,:), POINTER ::   U_zu          ! relative wind at zu                            [m/s]
      REAL(wp), DIMENSION(:,:), POINTER ::   Ce_n10        ! 10m neutral latent coefficient
      REAL(wp), DIMENSION(:,:), POINTER ::   Ch_n10        ! 10m neutral sensible coefficient
      REAL(wp), DIMENSION(:,:), POINTER ::   sqrt_Cd_n10   ! root square of Cd_n10
      REAL(wp), DIMENSION(:,:), POINTER ::   sqrt_Cd       ! root square of Cd
      REAL(wp), DIMENSION(:,:), POINTER ::   zeta_u        ! stability parameter at height zu
      REAL(wp), DIMENSION(:,:), POINTER ::   zeta_t        ! stability parameter at height zt
      REAL(wp), DIMENSION(:,:), POINTER ::   zpsi_h_u, zpsi_m_u
      REAL(wp), DIMENSION(:,:), POINTER ::   ztmp0, ztmp1, ztmp2
      REAL(wp), DIMENSION(:,:), POINTER ::   stab          ! 1st stability test integer
      !!----------------------------------------------------------------------

      IF( nn_timing == 1 )  CALL timing_start('turb_core_2z')
    
      CALL wrk_alloc( jpi,jpj, U_zu, Ce_n10, Ch_n10, sqrt_Cd_n10, sqrt_Cd )
      CALL wrk_alloc( jpi,jpj, zeta_u, stab )
      CALL wrk_alloc( jpi,jpj, zpsi_h_u, zpsi_m_u, ztmp0, ztmp1, ztmp2 )

      l_zt_equal_zu = .FALSE.
      IF( ABS(zu - zt) < 0.01 ) l_zt_equal_zu = .TRUE.    ! testing "zu == zt" is risky with double precision

      IF( .NOT. l_zt_equal_zu )   CALL wrk_alloc( jpi,jpj, zeta_t )

      U_zu = MAX( 0.5 , dU )   !  relative wind speed at zu (normally 10m), we don't want to fall under 0.5 m/s

      !! First guess of stability: 
      ztmp0 = T_zt*(1. + 0.608*q_zt) - sst*(1. + 0.608*q_sat) ! air-sea difference of virtual pot. temp. at zt
      stab  = 0.5 + sign(0.5,ztmp0)                           ! stab = 1 if dTv > 0  => STABLE, 0 if unstable

      !! Neutral coefficients at 10m:
      IF( ln_cdgw ) THEN      ! wave drag case
         cdn_wave(:,:) = cdn_wave(:,:) + rsmall * ( 1._wp - tmask(:,:,1) )
         ztmp0   (:,:) = cdn_wave(:,:)
      ELSE
         ztmp0 = cd_neutral_10m( U_zu )
      ENDIF
      sqrt_Cd_n10 = SQRT( ztmp0 )
      Ce_n10  = 1.e-3*( 34.6 * sqrt_Cd_n10 )
      Ch_n10  = 1.e-3*sqrt_Cd_n10*(18.*stab + 32.7*(1. - stab))
    
      !! Initializing transf. coeff. with their first guess neutral equivalents :
      Cd = ztmp0   ;   Ce = Ce_n10   ;   Ch = Ch_n10   ;   sqrt_Cd = sqrt_Cd_n10

      !! Initializing values at z_u with z_t values:
      T_zu = T_zt   ;   q_zu = q_zt

      !!  * Now starting iteration loop
      DO j_itt=1, nb_itt
         !
         ztmp1 = T_zu - sst   ! Updating air/sea differences
         ztmp2 = q_zu - q_sat 

         ! Updating turbulent scales :   (L&Y 2004 eq. (7))
         ztmp1  = Ch/sqrt_Cd*ztmp1    ! theta*
         ztmp2  = Ce/sqrt_Cd*ztmp2    ! q*
       
         ztmp0 = T_zu*(1. + 0.608*q_zu) ! virtual potential temperature at zu

         ! Estimate the inverse of Monin-Obukov length (1/L) at height zu:
         ztmp0 =  (vkarmn*grav/ztmp0*(ztmp1*(1.+0.608*q_zu) + 0.608*T_zu*ztmp2)) / (Cd*U_zu*U_zu) 
         !                                                                     ( Cd*U_zu*U_zu is U*^2 at zu)

         !! Stability parameters :
         zeta_u   = zu*ztmp0   ;  zeta_u = sign( min(abs(zeta_u),10.0), zeta_u )
         zpsi_h_u = psi_h( zeta_u )
         zpsi_m_u = psi_m( zeta_u )
       
         !! Shifting temperature and humidity at zu (L&Y 2004 eq. (9b-9c))
         IF ( .NOT. l_zt_equal_zu ) THEN
            zeta_t = zt*ztmp0 ;  zeta_t = sign( min(abs(zeta_t),10.0), zeta_t )
            stab = LOG(zu/zt) - zpsi_h_u + psi_h(zeta_t)  ! stab just used as temp array!!!
            T_zu = T_zt + ztmp1/vkarmn*stab    ! ztmp1 is still theta*
            q_zu = q_zt + ztmp2/vkarmn*stab    ! ztmp2 is still q*
            q_zu = max(0., q_zu)
         END IF
       
         IF( ln_cdgw ) THEN      ! surface wave case
            sqrt_Cd = vkarmn / ( vkarmn / sqrt_Cd_n10 - zpsi_m_u ) 
            Cd      = sqrt_Cd * sqrt_Cd
         ELSE
           ! Update neutral wind speed at 10m and neutral Cd at 10m (L&Y 2004 eq. 9a)...
           !   In very rare low-wind conditions, the old way of estimating the
           !   neutral wind speed at 10m leads to a negative value that causes the code
           !   to crash. To prevent this a threshold of 0.25m/s is imposed.
           ztmp0 = MAX( 0.25 , U_zu/(1. + sqrt_Cd_n10/vkarmn*(LOG(zu/10.) - zpsi_m_u)) ) !  U_n10
           ztmp0 = cd_neutral_10m(ztmp0)                                                 ! Cd_n10
           sqrt_Cd_n10 = sqrt(ztmp0)
       
           Ce_n10  = 1.e-3 * (34.6 * sqrt_Cd_n10)                     ! L&Y 2004 eq. (6b)
           stab    = 0.5 + sign(0.5,zeta_u)                           ! update stability
           Ch_n10  = 1.e-3*sqrt_Cd_n10*(18.*stab + 32.7*(1. - stab))  ! L&Y 2004 eq. (6c-6d)

           !! Update of transfer coefficients:
           ztmp1 = 1. + sqrt_Cd_n10/vkarmn*(LOG(zu/10.) - zpsi_m_u)   ! L&Y 2004 eq. (10a)
           Cd      = ztmp0 / ( ztmp1*ztmp1 )   
           sqrt_Cd = SQRT( Cd )
         ENDIF
         !
         ztmp0 = (LOG(zu/10.) - zpsi_h_u) / vkarmn / sqrt_Cd_n10
         ztmp2 = sqrt_Cd / sqrt_Cd_n10
         ztmp1 = 1. + Ch_n10*ztmp0               
         Ch  = Ch_n10*ztmp2 / ztmp1  ! L&Y 2004 eq. (10b)
         !
         ztmp1 = 1. + Ce_n10*ztmp0               
         Ce  = Ce_n10*ztmp2 / ztmp1  ! L&Y 2004 eq. (10c)
         !
      END DO

      CALL wrk_dealloc( jpi,jpj, U_zu, Ce_n10, Ch_n10, sqrt_Cd_n10, sqrt_Cd )
      CALL wrk_dealloc( jpi,jpj, zeta_u, stab )
      CALL wrk_dealloc( jpi,jpj, zpsi_h_u, zpsi_m_u, ztmp0, ztmp1, ztmp2 )

      IF( .NOT. l_zt_equal_zu ) CALL wrk_dealloc( jpi,jpj, zeta_t )

      IF( nn_timing == 1 )  CALL timing_stop('turb_core_2z')
      !
   END SUBROUTINE turb_core_2z


   FUNCTION cd_neutral_10m( zw10 )
      !!----------------------------------------------------------------------
      !! Estimate of the neutral drag coefficient at 10m as a function 
      !! of neutral wind  speed at 10m
      !!
      !! Origin: Large & Yeager 2008 eq.(11a) and eq.(11b)
      !!
      !! Author: L. Brodeau, june 2014
      !!----------------------------------------------------------------------    
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   zw10           ! scalar wind speed at 10m (m/s)
      REAL(wp), DIMENSION(jpi,jpj)             ::   cd_neutral_10m
      !
      REAL(wp), DIMENSION(:,:), POINTER ::   rgt33
      !!----------------------------------------------------------------------    
      !
      CALL wrk_alloc( jpi,jpj, rgt33 )
      !
      !! When wind speed > 33 m/s => Cyclone conditions => special treatment
      rgt33 = 0.5_wp + SIGN( 0.5_wp, (zw10 - 33._wp) )   ! If zw10 < 33. => 0, else => 1  
      cd_neutral_10m = 1.e-3 * ( &
         &       (1._wp - rgt33)*( 2.7_wp/zw10 + 0.142_wp + zw10/13.09_wp - 3.14807E-10*zw10**6) & ! zw10< 33.
         &      + rgt33         *      2.34   )                                                    ! zw10 >= 33.
      !
      CALL wrk_dealloc( jpi,jpj, rgt33)
      !
   END FUNCTION cd_neutral_10m


   FUNCTION psi_m(pta)   !! Psis, L&Y 2004 eq. (8c), (8d), (8e)
      !-------------------------------------------------------------------------------
      ! universal profile stability function for momentum
      !-------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pta
      !
      REAL(wp), DIMENSION(jpi,jpj)             :: psi_m
      REAL(wp), DIMENSION(:,:), POINTER        :: X2, X, stabit
      !-------------------------------------------------------------------------------
      !
      CALL wrk_alloc( jpi,jpj, X2, X, stabit )
      !
      X2 = SQRT( ABS( 1. - 16.*pta ) )  ;  X2 = MAX( X2 , 1. )   ;   X = SQRT( X2 )
      stabit = 0.5 + SIGN( 0.5 , pta )
      psi_m = -5.*pta*stabit  &                                                          ! Stable
         &    + (1. - stabit)*(2.*LOG((1. + X)*0.5) + LOG((1. + X2)*0.5) - 2.*ATAN(X) + rpi*0.5)  ! Unstable
      !
      CALL wrk_dealloc( jpi,jpj, X2, X, stabit )
      !
   END FUNCTION psi_m


   FUNCTION psi_h( pta )    !! Psis, L&Y 2004 eq. (8c), (8d), (8e)
      !-------------------------------------------------------------------------------
      ! universal profile stability function for temperature and humidity
      !-------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   pta
      !
      REAL(wp), DIMENSION(jpi,jpj)             ::   psi_h
      REAL(wp), DIMENSION(:,:), POINTER        ::   X2, X, stabit
      !-------------------------------------------------------------------------------
      !
      CALL wrk_alloc( jpi,jpj, X2, X, stabit )
      !
      X2 = SQRT( ABS( 1. - 16.*pta ) )   ;   X2 = MAX( X2 , 1. )   ;   X = SQRT( X2 )
      stabit = 0.5 + SIGN( 0.5 , pta )
      psi_h = -5.*pta*stabit   &                                       ! Stable
         &    + (1. - stabit)*(2.*LOG( (1. + X2)*0.5 ))                ! Unstable
      !
      CALL wrk_dealloc( jpi,jpj, X2, X, stabit )
      !
   END FUNCTION psi_h

   !!======================================================================
END MODULE sbcblk_core
