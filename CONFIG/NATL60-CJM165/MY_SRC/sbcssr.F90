MODULE sbcssr
   !!======================================================================
   !!                       ***  MODULE  sbcssr  ***
   !! Surface module :  heat and fresh water fluxes a restoring term toward observed SST/SSS
   !!======================================================================
   !! History :  3.0  !  2006-06  (G. Madec)  Original code
   !!            3.2  !  2009-04  (B. Lemaire)  Introduce iom_put
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   sbc_ssr       : add to sbc a restoring term toward SST/SSS climatology
   !!   sbc_ssr_init  : initialisation of surface restoring
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE sbc_oce        ! surface boundary condition
   USE phycst         ! physical constants
   USE sbcrnf         ! surface boundary condition : runoffs
   !
   USE fldread        ! read input fields
   USE iom            ! I/O manager
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! distribued memory computing library
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE timing         ! Timing
   USE shapiro         ! used in case of ln_sssr_flt
   USE lib_fortran    ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)  

   IMPLICIT NONE
   PRIVATE

   PUBLIC   sbc_ssr        ! routine called in sbcmod
   PUBLIC   sbc_ssr_init   ! routine called in sbcmod

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   erp   !: evaporation damping   [kg/m2/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   qrp   !: heat flux damping        [w/m2]
   !{ DRAKKAR : local increase of restoring. eg Med Sea in ORCA025
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   erpcoef  !: multiplicating coef for local change to erp
   !}

   !                                   !!* Namelist namsbc_ssr *
   INTEGER, PUBLIC ::   nn_sstr         ! SST/SSS restoring indicator
   INTEGER, PUBLIC ::   nn_sssr         ! SST/SSS restoring indicator
   REAL(wp)        ::   rn_dqdt         ! restoring factor on SST and SSS
   REAL(wp)        ::   rn_deds         ! restoring factor on SST and SSS
   LOGICAL         ::   ln_sssr_bnd     ! flag to bound erp term 
   REAL(wp)        ::   rn_sssr_bnd     ! ABS(Max./Min.) value of erp term [mm/day]
   !{ DRAKKAR 
   LOGICAL, PUBLIC ::   ln_sssr_flt     ! flag to filter sss for restoring
   INTEGER, PUBLIC ::   nn_shap_iter    ! number of iteration for shapiro
   !}

   REAL(wp) , ALLOCATABLE, DIMENSION(:) ::   buffer   ! Temporary buffer for exchange
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_sst   ! structure of input SST (file informations, fields read)
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_sss   ! structure of input SSS (file informations, fields read)
   !{ DRAKKAR : limit sss restoring in the coastal area
   LOGICAL         :: ln_sssr_msk
   TYPE(FLD_N)     :: sn_coast
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: distcoast   ! use to read the distance and then for weight purpose

   REAL(wp)        :: rn_dist      ! (km) decaying lenght scale for SSS restoring near the coast
   !}

   !! * Substitutions
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO Consortium (2011)
   !! $Id: sbcssr.F90 4990 2014-12-15 16:42:49Z timgraham $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE sbc_ssr( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE sbc_ssr  ***
      !!
      !! ** Purpose :   Add to heat and/or freshwater fluxes a damping term
      !!                toward observed SST and/or SSS.
      !!
      !! ** Method  : - Read namelist namsbc_ssr
      !!              - Read observed SST and/or SSS
      !!              - at each nscb time step
      !!                   add a retroaction term on qns    (nn_sstr = 1)
      !!                   add a damping term on sfx        (nn_sssr = 1)
      !!                   add a damping term on emp        (nn_sssr = 2)
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in   ) ::   kt   ! ocean time step
      !!
      INTEGER  ::   ji, jj   ! dummy loop indices
      REAL(wp) ::   zerp     ! local scalar for evaporation damping
      REAL(wp) ::   zqrp     ! local scalar for heat flux damping
      REAL(wp) ::   zsrp     ! local scalar for unit conversion of rn_deds factor
      REAL(wp) ::   zerp_bnd ! local scalar for unit conversion of rn_epr_max factor
      INTEGER  ::   ierror   ! return error code
      !{ DRAKKAR 
      REAL(wp) , DIMENSION (jpi,jpj) :: zsss_m ! temporary array
      !}
      !!
      CHARACTER(len=100) ::  cn_dir          ! Root directory for location of ssr files
      TYPE(FLD_N) ::   sn_sst, sn_sss        ! informations about the fields to be read
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('sbc_ssr')
      !
      IF( nn_sstr + nn_sssr /= 0 ) THEN
         !
         IF( nn_sstr == 1)   CALL fld_read( kt, nn_fsbc, sf_sst )   ! Read SST data and provides it at kt
         IF( nn_sssr >= 1)   CALL fld_read( kt, nn_fsbc, sf_sss )   ! Read SSS data and provides it at kt
         !
         !                                         ! ========================= !
         IF( MOD( kt-1, nn_fsbc ) == 0 ) THEN      !    Add restoring term     !
            !                                      ! ========================= !
            !
            IF( nn_sstr == 1 ) THEN                                   !* Temperature restoring term
               DO jj = 1, jpj
                  DO ji = 1, jpi
                     zqrp = rn_dqdt * ( sst_m(ji,jj) - sf_sst(1)%fnow(ji,jj,1) )
                     qns(ji,jj) = qns(ji,jj) + zqrp
                     qrp(ji,jj) = zqrp
                  END DO
               END DO
               CALL iom_put( "qrp", qrp )                             ! heat flux damping
            ENDIF
            !
            IF( nn_sssr == 1 ) THEN                                   !* Salinity damping term (salt flux only (sfx))
               zsrp = rn_deds / rday                                  ! from [mm/day] to [kg/m2/s]
!CDIR COLLAPSE
               DO jj = 1, jpj
                  DO ji = 1, jpi
                     zerp = zsrp * ( 1. - 2.*rnfmsk(ji,jj) )   &      ! No damping in vicinity of river mouths
                        &        * ( sss_m(ji,jj) - sf_sss(1)%fnow(ji,jj,1) ) * erpcoef(ji,jj)
                     sfx(ji,jj) = sfx(ji,jj) + zerp                 ! salt flux
                     erp(ji,jj) = zerp / MAX( sss_m(ji,jj), 1.e-20 ) ! converted into an equivalent volume flux (diagnostic only)
                  END DO
               END DO
               CALL iom_put( "erp", erp )                             ! freshwater flux damping
               !
            ELSEIF( nn_sssr == 2 ) THEN                               !* Salinity damping term (volume flux (emp) and associated heat flux (qns)
               zsrp = rn_deds / rday                                  ! from [mm/day] to [kg/m2/s]
               zerp_bnd = rn_sssr_bnd / rday                          !       -              -    
               !{ DRAKKAR : using filtered sss for restoring 
               IF (ln_sssr_flt ) THEN
                  CALL Shapiro_1D ( sss_m(:,:), nn_shap_iter, 'ORCA_GLOB', zsss_m )
                  zsss_m = zsss_m * tmask(:,:,1)
               ELSE
                  zsss_m = sss_m * tmask(:,:,1)
               ENDIF
               !}
!CDIR COLLAPSE
               DO jj = 1, jpj
                  DO ji = 1, jpi                            
                     zerp = zsrp * ( 1. - 2.*rnfmsk(ji,jj) )   &      ! No damping in vicinity of river mouths
                        &        * ( sss_m(ji,jj) - sf_sss(1)%fnow(ji,jj,1) )   &
                        &        / MAX(  sss_m(ji,jj), 1.e-20   )               &
                          &        * erpcoef(ji,jj)
                     IF( ln_sssr_bnd )   zerp = SIGN( 1., zerp ) * MIN( zerp_bnd, ABS(zerp) )
                     IF( ln_sssr_msk )   zerp = zerp * distcoast(ji,jj) ! multiply by weigh to fade zerp out near the coast
                     emp (ji,jj) = emp (ji,jj) + zerp
                     qns(ji,jj) = qns(ji,jj) - zerp * rcp * sst_m(ji,jj)
                     erp (ji,jj) = zerp
                  END DO
               END DO
               CALL iom_put( "erp", erp )                             ! freshwater flux damping
            ENDIF
            !
         ENDIF
         !
      ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop('sbc_ssr')
      !
   END SUBROUTINE sbc_ssr

 
   SUBROUTINE sbc_ssr_init
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE sbc_ssr_init  ***
      !!
      !! ** Purpose :   initialisation of surface damping term
      !!
      !! ** Method  : - Read namelist namsbc_ssr
      !!              - Read observed SST and/or SSS if required
      !!---------------------------------------------------------------------
      INTEGER  ::   ji, jj   ! dummy loop indices
      REAL(wp) ::   zerp     ! local scalar for evaporation damping
      REAL(wp) ::   zqrp     ! local scalar for heat flux damping
      REAL(wp) ::   zsrp     ! local scalar for unit conversion of rn_deds factor
      REAL(wp) ::   zerp_bnd ! local scalar for unit conversion of rn_epr_max factor
      INTEGER  ::   ierror   ! return error code
      !!
      !{ DRAKKAR 
      INTEGER  ::   ii0, ii1, ii2, ij0, ij1, ij2, inum
      REAL(wp) :: zalph
      !}
      CHARACTER(len=100) ::  cn_dir          ! Root directory for location of ssr files
      CHARACTER(LEN=100) ::  cl_coastfile
      TYPE(FLD_N) ::   sn_sst, sn_sss        ! informations about the fields to be read
      NAMELIST/namsbc_ssr/ cn_dir, nn_sstr, nn_sssr, rn_dqdt, rn_deds, sn_sst, sn_sss
      NAMELIST/namsbc_ssr/ ln_sssr_bnd, rn_sssr_bnd, ln_sssr_flt, ln_sssr_msk, sn_coast, rn_dist, nn_shap_iter
      INTEGER     ::  ios
      !!----------------------------------------------------------------------
      !
 
      REWIND( numnam_ref )              ! Namelist namsbc_ssr in reference namelist : 
      READ  ( numnam_ref, namsbc_ssr, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namsbc_ssr in reference namelist', lwp )

      REWIND( numnam_cfg )              ! Namelist namsbc_ssr in configuration namelist :
      READ  ( numnam_cfg, namsbc_ssr, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namsbc_ssr in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, namsbc_ssr )

      IF(lwp) THEN                 !* control print
         WRITE(numout,*)
         WRITE(numout,*) 'sbc_ssr : SST and/or SSS damping term '
         WRITE(numout,*) '~~~~~~~ '
         WRITE(numout,*) '   Namelist namsbc_ssr :'
         WRITE(numout,*) '      SST restoring term (Yes=1)             nn_sstr     = ', nn_sstr
         WRITE(numout,*) '      SSS damping term (Yes=1, salt flux)    nn_sssr     = ', nn_sssr
         WRITE(numout,*) '                       (Yes=2, volume flux) '
         WRITE(numout,*) '      dQ/dT (restoring magnitude on SST)     rn_dqdt     = ', rn_dqdt, ' W/m2/K'
         WRITE(numout,*) '      dE/dS (restoring magnitude on SST)     rn_deds     = ', rn_deds, ' mm/day'
         WRITE(numout,*) '      flag to bound erp term                 ln_sssr_bnd = ', ln_sssr_bnd
         WRITE(numout,*) '      ABS(Max./Min.) erp threshold           rn_sssr_bnd = ', rn_sssr_bnd, ' mm/day'
         WRITE(numout,*) '      Filtering of sss for restoring         ln_sssr_flt = ', ln_sssr_flt 
         IF ( ln_sssr_flt ) THEN
            WRITE(numout,*) '      Number of used Shapiro filter           nn_shap_iter = ', nn_shap_iter
         ENDIF
         WRITE(numout,*) '      Limit sss restoring near the coast     ln_sssr_msk = ', ln_sssr_msk
         IF ( ln_sssr_msk ) WRITE(numout,*) '      Decaying lenght scale from the coast   rn_dist     = ', rn_dist, ' km'
      ENDIF
      !
      !                            !* Allocate erp and qrp array
      ALLOCATE( qrp(jpi,jpj), erp(jpi,jpj), erpcoef(jpi,jpj),STAT=ierror )
      IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_ssr: unable to allocate erp and qrp array' )
      !
      IF( nn_sstr == 1 ) THEN      !* set sf_sst structure & allocate arrays
         !
         ALLOCATE( sf_sst(1), STAT=ierror )
         IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_ssr: unable to allocate sf_sst structure' )
         ALLOCATE( sf_sst(1)%fnow(jpi,jpj,1), STAT=ierror )
         IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_ssr: unable to allocate sf_sst now array' )
         !
         ! fill sf_sst with sn_sst and control print
         CALL fld_fill( sf_sst, (/ sn_sst /), cn_dir, 'sbc_ssr', 'SST restoring term toward SST data', 'namsbc_ssr' )
         IF( sf_sst(1)%ln_tint )   ALLOCATE( sf_sst(1)%fdta(jpi,jpj,1,2), STAT=ierror )
         IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_ssr: unable to allocate sf_sst data array' )
         !
      ENDIF
      !
      IF( nn_sssr >= 1 ) THEN      !* set sf_sss structure & allocate arrays
         !
         ALLOCATE( sf_sss(1), STAT=ierror )
         IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_ssr: unable to allocate sf_sss structure' )
         ALLOCATE( sf_sss(1)%fnow(jpi,jpj,1), STAT=ierror )
         IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_ssr: unable to allocate sf_sss now array' )
         !
         ! fill sf_sss with sn_sss and control print
         CALL fld_fill( sf_sss, (/ sn_sss /), cn_dir, 'sbc_ssr', 'SSS restoring term toward SSS data', 'namsbc_ssr' )
         IF( sf_sss(1)%ln_tint )   ALLOCATE( sf_sss(1)%fdta(jpi,jpj,1,2), STAT=ierror )
         IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_ssr: unable to allocate sf_sss data array' )
         !
         ! if masking of coastal area is used
         IF ( ln_sssr_msk ) THEN
            ALLOCATE( distcoast(jpi,jpj),STAT=ierror )  
            IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_ssr: unable to allocate erp and qrp array' )
            WRITE(cl_coastfile,'(a,a)' ) TRIM( cn_dir ), TRIM( sn_coast%clname )
            CALL iom_open ( cl_coastfile, inum )                          ! open file
            CALL iom_get  ( inum, jpdom_data, sn_coast%clvar, distcoast ) ! read tcoast  in m
            CALL iom_close( inum )
            ! transform distcoast to weight 
            rn_dist=rn_dist*1000.  ! tranform rn_dist to m
            distcoast(:,:)=0.5*(tanh(3.*(distcoast(:,:)*distcoast(:,:)/rn_dist/rn_dist - 1 )) + 1 )
         ENDIF
      ENDIF
      !
      ! DRAKKAR { initialize erpcoef to increase erp in the med sea
      erpcoef(:,:) = 1._wp
      IF( cp_cfg == "orca" .AND. jp_cfg == 25 ) THEN  ! ORCA R025 configuration
         !! add extra SSS restoring in the Red Sea
         ii0= 1280       ;  ii1 = 1325
         ij0= 560        ;  ij1 = 625
         erpcoef( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1)) =  5.0

         !! add extra SSS restoring in the Med sea (x3) decreasing to 1 into the alboran sea
         ii1 = 1145  ; ii2 =  1330
         ij1 = 626   ; ij2 =  726
         erpcoef( mi0(ii1):mi1(ii2) , mj0(ij1):mj1(ij2)) = 3.0

         !! decrease in alboran sea (along i )
         ii0= 1128       ;  ii1 = 1144
         ij0= 645        ;  ij1 = 670
         DO jj=mj0(ij0), mj1(ij1)
            DO ji= mi0(ii0), mi1(ii1)
               !zalph=( alph1 -alph0 )* (I - ii0 )/(ii1-ii0) + alph0
               !zalph=( 3.    - 1.   )* (I - ii0 )/(ii1-ii0) + 1.
               zalph=     2.   *    (mig(ji)-ii0)/(ii1-ii0) + 1.
               erpcoef(ji,jj) = zalph 
            ENDDO
         ENDDO
      ENDIF
      !            }
      !                            !* Initialize qrp and erp if no restoring 
      IF( nn_sstr /= 1                   )   qrp(:,:) = 0._wp
      IF( nn_sssr /= 1 .OR. nn_sssr /= 2 )   erp(:,:) = 0._wp
      !
   END SUBROUTINE sbc_ssr_init
      
   !!======================================================================
END MODULE sbcssr
