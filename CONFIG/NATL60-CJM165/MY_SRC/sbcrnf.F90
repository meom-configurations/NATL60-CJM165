MODULE sbcrnf
   !!======================================================================
   !!                       ***  MODULE  sbcrnf  ***
   !! Ocean forcing:  river runoff
   !!=====================================================================
   !! History :  OPA  ! 2000-11  (R. Hordoir, E. Durand)  NetCDF FORMAT
   !!   NEMO     1.0  ! 2002-09  (G. Madec)  F90: Free form and module
   !!            3.0  ! 2006-07  (G. Madec)  Surface module
   !!            3.2  ! 2009-04  (B. Lemaire)  Introduce iom_put
   !!            3.3  ! 2010-10  (R. Furner, G. Madec) runoff distributed over ocean levels
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   sbc_rnf      : monthly runoffs read in a NetCDF file
   !!   sbc_rnf_init : runoffs initialisation
   !!   rnf_mouth    : set river mouth mask
   !!----------------------------------------------------------------------
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE sbc_oce         ! surface boundary condition variables
   USE sbcisf          ! PM we could remove it I think
   USE closea          ! closed seas
   USE fldread         ! read input field at current time step
   USE in_out_manager  ! I/O manager
   USE iom             ! I/O module
   USE lib_mpp         ! MPP library
   USE eosbn2
   USE wrk_nemo        ! Memory allocation

   IMPLICIT NONE
   PRIVATE

   PUBLIC   sbc_rnf       ! routine call in sbcmod module
   PUBLIC   sbc_rnf_div   ! routine called in divcurl module
   PUBLIC   sbc_rnf_alloc ! routine call in sbcmod module
   PUBLIC   sbc_rnf_init  ! (PUBLIC for TAM)
   !                                                     !!* namsbc_rnf namelist *
   CHARACTER(len=100)         ::   cn_dir          !: Root directory for location of rnf files
   LOGICAL                    ::   ln_rnf_depth      !: depth       river runoffs attribute specified in a file
   LOGICAL                    ::   ln_rnf_depth_ini  !: depth       river runoffs  computed at the initialisation
   REAL(wp)                   ::   rn_rnf_max        !: maximum value of the runoff climatologie ( ln_rnf_depth_ini = .true )
   REAL(wp)                   ::   rn_dep_max        !: depth over which runoffs is spread ( ln_rnf_depth_ini = .true )
   INTEGER                    ::   nn_rnf_depth_file !: create (=1) a runoff depth file or not (=0)
   LOGICAL                    ::   ln_rnf_tem      !: temperature river runoffs attribute specified in a file
   LOGICAL           , PUBLIC ::   ln_rnf_sal      !: salinity    river runoffs attribute specified in a file
   TYPE(FLD_N)       , PUBLIC ::   sn_rnf          !: information about the runoff file to be read
   TYPE(FLD_N)                ::   sn_cnf          !: information about the runoff mouth file to be read
   TYPE(FLD_N)                ::   sn_s_rnf        !: information about the salinities of runoff file to be read
   TYPE(FLD_N)                ::   sn_t_rnf        !: information about the temperatures of runoff file to be read
   TYPE(FLD_N)                ::   sn_dep_rnf      !: information about the depth which river inflow affects
   LOGICAL           , PUBLIC ::   ln_rnf_mouth    !: specific treatment in mouths vicinity
   REAL(wp)                   ::   rn_hrnf         !: runoffs, depth over which enhanced vertical mixing is used
   REAL(wp)          , PUBLIC ::   rn_avt_rnf      !: runoffs, value of the additional vertical mixing coef. [m2/s]
   REAL(wp)          , PUBLIC ::   rn_rfact        !: multiplicative factor for runoff

   LOGICAL           , PUBLIC ::   l_rnfcpl = .false.       ! runoffs recieved from oasis

   INTEGER , PUBLIC  ::   nkrnf = 0         !: nb of levels over which Kz is increased at river mouths
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   rnfmsk              !: river mouth mask (hori.)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)     ::   rnfmsk_z            !: river mouth mask (vert.)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   h_rnf               !: depth of runoff in m
   INTEGER,  PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   nk_rnf              !: depth of runoff in model levels
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   rnf_tsc_b, rnf_tsc  !: before and now T & S runoff contents   [K.m/s & PSU.m/s]   

   TYPE(FLD),        ALLOCATABLE, DIMENSION(:) ::   sf_rnf       ! structure: river runoff (file information, fields read)
   TYPE(FLD),        ALLOCATABLE, DIMENSION(:) ::   sf_s_rnf     ! structure: river runoff salinity (file information, fields read)  
   TYPE(FLD),        ALLOCATABLE, DIMENSION(:) ::   sf_t_rnf     ! structure: river runoff temperature (file information, fields read)  
 
   !! * Substitutions  
#  include "domzgr_substitute.h90"  
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: sbcrnf.F90 5503 2015-06-29 12:31:29Z cetlod $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION sbc_rnf_alloc()
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE sbc_rnf_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( rnfmsk(jpi,jpj)         , rnfmsk_z(jpk)          ,     &
         &      h_rnf (jpi,jpj)         , nk_rnf  (jpi,jpj)      ,     &
         &      rnf_tsc_b(jpi,jpj,jpts) , rnf_tsc (jpi,jpj,jpts) , STAT=sbc_rnf_alloc )
         !
      IF( lk_mpp            )   CALL mpp_sum ( sbc_rnf_alloc )
      IF( sbc_rnf_alloc > 0 )   CALL ctl_warn('sbc_rnf_alloc: allocation of arrays failed')
   END FUNCTION sbc_rnf_alloc


   SUBROUTINE sbc_rnf( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sbc_rnf  ***
      !!
      !! ** Purpose :   Introduce a climatological run off forcing
      !!
      !! ** Method  :   Set each river mouth with a monthly climatology
      !!                provided from different data.
      !!                CAUTION : upward water flux, runoff forced to be < 0
      !!
      !! ** Action  :   runoff updated runoff field at time-step kt
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt          ! ocean time step
      !
      INTEGER  ::   ji, jj    ! dummy loop indices
      INTEGER  ::   z_err = 0 ! dummy integer for error handling
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), POINTER       ::   ztfrz   ! freezing point used for temperature correction
      !
      CALL wrk_alloc( jpi,jpj, ztfrz)

      !                                            ! ---------------------------------------- !
      IF( kt /= nit000 ) THEN                      !          Swap of forcing fields          !
         !                                         ! ---------------------------------------- !
         rnf_b    (:,:  ) = rnf    (:,:  )               ! Swap the ocean forcing fields except at nit000
         rnf_tsc_b(:,:,:) = rnf_tsc(:,:,:)               ! where before fields are set at the end of the routine
         !
      ENDIF

      !                                            !-------------------!
      !                                            !   Update runoff   !
      !                                            !-------------------!
      !
      IF( .NOT. l_rnfcpl )   CALL fld_read ( kt, nn_fsbc, sf_rnf   )    ! Read Runoffs data and provide it at kt
      IF(   ln_rnf_tem   )   CALL fld_read ( kt, nn_fsbc, sf_t_rnf )    ! idem for runoffs temperature if required
      IF(   ln_rnf_sal   )   CALL fld_read ( kt, nn_fsbc, sf_s_rnf )    ! idem for runoffs salinity    if required
      !
      IF( MOD( kt - 1, nn_fsbc ) == 0 ) THEN
         !
         IF( .NOT. l_rnfcpl )   rnf(:,:) = rn_rfact * ( sf_rnf(1)%fnow(:,:,1) )       ! updated runoff value at time step kt
         !
         !                                                     ! set temperature & salinity content of runoffs
         IF( ln_rnf_tem ) THEN                                       ! use runoffs temperature data
            rnf_tsc(:,:,jp_tem) = ( sf_t_rnf(1)%fnow(:,:,1) ) * rnf(:,:) * r1_rau0
            WHERE( sf_t_rnf(1)%fnow(:,:,1) == -999._wp )             ! if missing data value use SST as runoffs temperature
               rnf_tsc(:,:,jp_tem) = sst_m(:,:) * rnf(:,:) * r1_rau0
            END WHERE
            WHERE( sf_t_rnf(1)%fnow(:,:,1) == -222._wp )             ! where fwf comes from melting of ice shelves or iceberg
               ztfrz(:,:) = -1.9 !tfreez( sss_m(:,:) ) !PM to be discuss (trouble if sensitivity study)
               rnf_tsc(:,:,jp_tem) = ztfrz(:,:) * rnf(:,:) * r1_rau0 - rnf(:,:) * lfusisf * r1_rau0_rcp
            END WHERE
         ELSE                                                        ! use SST as runoffs temperature
            rnf_tsc(:,:,jp_tem) = sst_m(:,:) * rnf(:,:) * r1_rau0
         ENDIF
         !                                                           ! use runoffs salinity data
         IF( ln_rnf_sal )   rnf_tsc(:,:,jp_sal) = ( sf_s_rnf(1)%fnow(:,:,1) ) * rnf(:,:) * r1_rau0
         !                                                           ! else use S=0 for runoffs (done one for all in the init)
         CALL iom_put( "runoffs", rnf )         ! output runoffs arrays
      ENDIF
      !
      !                                                ! ---------------------------------------- !
      IF( kt == nit000 ) THEN                          !   set the forcing field at nit000 - 1    !
         !                                             ! ---------------------------------------- !
         IF( ln_rstart .AND.    &                               !* Restart: read in restart file
            & iom_varid( numror, 'rnf_b', ldstop = .FALSE. ) > 0 ) THEN
            IF(lwp) WRITE(numout,*) '          nit000-1 runoff forcing fields red in the restart file'
            CALL iom_get( numror, jpdom_autoglo, 'rnf_b', rnf_b )     ! before runoff
            CALL iom_get( numror, jpdom_autoglo, 'rnf_hc_b', rnf_tsc_b(:,:,jp_tem) )   ! before heat content of runoff
            CALL iom_get( numror, jpdom_autoglo, 'rnf_sc_b', rnf_tsc_b(:,:,jp_sal) )   ! before salinity content of runoff
         ELSE                                                   !* no restart: set from nit000 values
            IF(lwp) WRITE(numout,*) '          nit000-1 runoff forcing fields set to nit000'
            rnf_b    (:,:  ) = rnf    (:,:  )
            rnf_tsc_b(:,:,:) = rnf_tsc(:,:,:)
         ENDIF
      ENDIF
      !                                                ! ---------------------------------------- !
      IF( lrst_oce ) THEN                              !      Write in the ocean restart file     !
         !                                             ! ---------------------------------------- !
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'sbcrnf : runoff forcing fields written in ocean restart file ',   &
            &                    'at it= ', kt,' date= ', ndastp
         IF(lwp) WRITE(numout,*) '~~~~'
         CALL iom_rstput( kt, nitrst, numrow, 'rnf_b' , rnf )
         CALL iom_rstput( kt, nitrst, numrow, 'rnf_hc_b', rnf_tsc(:,:,jp_tem) )
         CALL iom_rstput( kt, nitrst, numrow, 'rnf_sc_b', rnf_tsc(:,:,jp_sal) )
      ENDIF
      !
      CALL wrk_dealloc( jpi,jpj, ztfrz)
      !
   END SUBROUTINE sbc_rnf


   SUBROUTINE sbc_rnf_div( phdivn )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sbc_rnf  ***
      !!
      !! ** Purpose :   update the horizontal divergence with the runoff inflow
      !!
      !! ** Method  :
      !!                CAUTION : rnf is positive (inflow) decreasing the
      !!                          divergence and expressed in m/s
      !!
      !! ** Action  :   phdivn   decreased by the runoff inflow
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) ::   phdivn   ! horizontal divergence
      !!
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   zfact     ! local scalar
      !!----------------------------------------------------------------------
      !
      zfact = 0.5_wp
      !
      IF( ln_rnf_depth .OR. ln_rnf_depth_ini ) THEN      !==   runoff distributed over several levels   ==!
         IF( lk_vvl ) THEN             ! variable volume case
            DO jj = 1, jpj                   ! update the depth over which runoffs are distributed
               DO ji = 1, jpi
                  h_rnf(ji,jj) = 0._wp
                  DO jk = 1, nk_rnf(ji,jj)                           ! recalculates h_rnf to be the depth in metres
                     h_rnf(ji,jj) = h_rnf(ji,jj) + fse3t(ji,jj,jk)   ! to the bottom of the relevant grid box
                  END DO
                  !                          ! apply the runoff input flow
                  DO jk = 1, nk_rnf(ji,jj)
                     phdivn(ji,jj,jk) = phdivn(ji,jj,jk) - ( rnf(ji,jj) + rnf_b(ji,jj) ) * zfact * r1_rau0 / h_rnf(ji,jj)
                  END DO
               END DO
            END DO
         ELSE                          ! constant volume case : just apply the runoff input flow
            DO jj = 1, jpj
               DO ji = 1, jpi
                  DO jk = 1, nk_rnf(ji,jj)
                     phdivn(ji,jj,jk) = phdivn(ji,jj,jk) - ( rnf(ji,jj) + rnf_b(ji,jj) ) * zfact * r1_rau0 / h_rnf(ji,jj)
                  END DO
               END DO
            END DO
         ENDIF
      ELSE                       !==   runoff put only at the surface   ==!
         IF( lk_vvl ) THEN              ! variable volume case
            h_rnf(:,:) = fse3t(:,:,1)   ! recalculate h_rnf to be depth of top box
         ENDIF
         phdivn(:,:,1) = phdivn(:,:,1) - ( rnf(:,:) + rnf_b(:,:) ) * zfact * r1_rau0 / fse3t(:,:,1)
      ENDIF
      !
   END SUBROUTINE sbc_rnf_div


   SUBROUTINE sbc_rnf_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sbc_rnf_init  ***
      !!
      !! ** Purpose :   Initialisation of the runoffs if (ln_rnf=T)
      !!
      !! ** Method  : - read the runoff namsbc_rnf namelist
      !!
      !! ** Action  : - read parameters
      !!----------------------------------------------------------------------
      CHARACTER(len=32) ::   rn_dep_file   ! runoff file name
      INTEGER           ::   ji, jj, jk, jm    ! dummy loop indices
      INTEGER           ::   ierror, inum  ! temporary integer
      INTEGER           ::   ios           ! Local integer output status for namelist read
      INTEGER           ::   nbrec         ! temporary integer
      REAL(wp)          ::   zacoef  
      REAL(wp), DIMENSION(12)                 :: zrec             ! times records
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: zrnfcl    
      REAL(wp), DIMENSION(:,:  ), ALLOCATABLE :: zrnf
      !
      NAMELIST/namsbc_rnf/ cn_dir            , ln_rnf_depth, ln_rnf_tem, ln_rnf_sal,   &
         &                 sn_rnf, sn_cnf    , sn_s_rnf    , sn_t_rnf  , sn_dep_rnf,   &
         &                 ln_rnf_mouth      , rn_hrnf     , rn_avt_rnf, rn_rfact,     &
         &                 ln_rnf_depth_ini  , rn_dep_max  , rn_rnf_max, nn_rnf_depth_file
      !!----------------------------------------------------------------------
      !
      !                                         !==  allocate runoff arrays
      IF( sbc_rnf_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'sbc_rnf_alloc : unable to allocate arrays' )
      !
      IF( .NOT. ln_rnf ) THEN                      ! no specific treatment in vicinity of river mouths 
         ln_rnf_mouth  = .FALSE.                   ! default definition needed for example by sbc_ssr or by tra_adv_muscl
         nkrnf         = 0
         rnf     (:,:) = 0.0_wp
         rnf_b   (:,:) = 0.0_wp
         rnfmsk  (:,:) = 0.0_wp
         rnfmsk_z(:)   = 0.0_wp
         RETURN
      ENDIF
      !
      !                                   ! ============
      !                                   !   Namelist
      !                                   ! ============
      !
      REWIND( numnam_ref )              ! Namelist namsbc_rnf in reference namelist : Runoffs 
      READ  ( numnam_ref, namsbc_rnf, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namsbc_rnf in reference namelist', lwp )

      REWIND( numnam_cfg )              ! Namelist namsbc_rnf in configuration namelist : Runoffs
      READ  ( numnam_cfg, namsbc_rnf, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namsbc_rnf in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, namsbc_rnf )
      !
      !                                         ! Control print
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'sbc_rnf : runoff '
         WRITE(numout,*) '~~~~~~~ '
         WRITE(numout,*) '   Namelist namsbc_rnf'
         WRITE(numout,*) '      specific river mouths treatment            ln_rnf_mouth = ', ln_rnf_mouth
         WRITE(numout,*) '      river mouth additional Kz                  rn_avt_rnf   = ', rn_avt_rnf
         WRITE(numout,*) '      depth of river mouth additional mixing     rn_hrnf      = ', rn_hrnf
         WRITE(numout,*) '      multiplicative factor for runoff           rn_rfact     = ', rn_rfact
      ENDIF
      !                                   ! ==================
      !                                   !   Type of runoff
      !                                   ! ==================
      !
      IF( .NOT. l_rnfcpl ) THEN                    
         ALLOCATE( sf_rnf(1), STAT=ierror )         ! Create sf_rnf structure (runoff inflow)
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '          runoffs inflow read in a file'
         IF( ierror > 0 ) THEN
            CALL ctl_stop( 'sbc_rnf: unable to allocate sf_rnf structure' )   ;   RETURN
         ENDIF
         ALLOCATE( sf_rnf(1)%fnow(jpi,jpj,1)   )
         IF( sn_rnf%ln_tint ) ALLOCATE( sf_rnf(1)%fdta(jpi,jpj,1,2) )
         CALL fld_fill( sf_rnf, (/ sn_rnf /), cn_dir, 'sbc_rnf_init', 'read runoffs data', 'namsbc_rnf' )
      ENDIF
      !
      IF( ln_rnf_tem ) THEN                      ! Create (if required) sf_t_rnf structure
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '          runoffs temperatures read in a file'
         ALLOCATE( sf_t_rnf(1), STAT=ierror  )
         IF( ierror > 0 ) THEN
            CALL ctl_stop( 'sbc_rnf_init: unable to allocate sf_t_rnf structure' )   ;   RETURN
         ENDIF
         ALLOCATE( sf_t_rnf(1)%fnow(jpi,jpj,1)   )
         IF( sn_t_rnf%ln_tint ) ALLOCATE( sf_t_rnf(1)%fdta(jpi,jpj,1,2) )
         CALL fld_fill (sf_t_rnf, (/ sn_t_rnf /), cn_dir, 'sbc_rnf_init', 'read runoff temperature data', 'namsbc_rnf' )
      ENDIF
      !
      IF( ln_rnf_sal  ) THEN                     ! Create (if required) sf_s_rnf and sf_t_rnf structures
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '          runoffs salinities read in a file'
         ALLOCATE( sf_s_rnf(1), STAT=ierror  )
         IF( ierror > 0 ) THEN
            CALL ctl_stop( 'sbc_rnf_init: unable to allocate sf_s_rnf structure' )   ;   RETURN
         ENDIF
         ALLOCATE( sf_s_rnf(1)%fnow(jpi,jpj,1)   )
         IF( sn_s_rnf%ln_tint ) ALLOCATE( sf_s_rnf(1)%fdta(jpi,jpj,1,2) )
         CALL fld_fill (sf_s_rnf, (/ sn_s_rnf /), cn_dir, 'sbc_rnf_init', 'read runoff salinity data', 'namsbc_rnf' )
      ENDIF
      !
      IF( ln_rnf_depth ) THEN                    ! depth of runoffs set from a file
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '          runoffs depth read in a file'
         rn_dep_file = TRIM( cn_dir )//TRIM( sn_dep_rnf%clname )
         IF( .NOT. sn_dep_rnf%ln_clim ) THEN   ;   WRITE(rn_dep_file, '(a,"_y",i4)' ) TRIM( rn_dep_file ), nyear    ! add year 
            IF( sn_dep_rnf%cltype == 'monthly' )   WRITE(rn_dep_file, '(a,"m",i2)'  ) TRIM( rn_dep_file ), nmonth   ! add month 
         ENDIF
         CALL iom_open ( rn_dep_file, inum )                           ! open file
         CALL iom_get  ( inum, jpdom_data, sn_dep_rnf%clvar, h_rnf )   ! read the river mouth array
         CALL iom_close( inum )                                        ! close file
         !
         nk_rnf(:,:) = 0                               ! set the number of level over which river runoffs are applied
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( h_rnf(ji,jj) > 0._wp ) THEN
                  jk = 2
                  DO WHILE ( jk /= mbkt(ji,jj) .AND. gdept_0(ji,jj,jk) < h_rnf(ji,jj) ) ;  jk = jk + 1
                  END DO
                  nk_rnf(ji,jj) = jk
               ELSEIF( h_rnf(ji,jj) == -1._wp   ) THEN   ;  nk_rnf(ji,jj) = 1
               ELSEIF( h_rnf(ji,jj) == -999._wp ) THEN   ;  nk_rnf(ji,jj) = mbkt(ji,jj)
               ELSE
                  CALL ctl_stop( 'sbc_rnf_init: runoff depth not positive, and not -999 or -1, rnf value in file fort.999'  )
                  WRITE(999,*) 'ji, jj, h_rnf(ji,jj) :', ji, jj, h_rnf(ji,jj)
               ENDIF
            END DO
         END DO
         DO jj = 1, jpj                                ! set the associated depth
            DO ji = 1, jpi
               h_rnf(ji,jj) = 0._wp
               DO jk = 1, nk_rnf(ji,jj)
                  h_rnf(ji,jj) = h_rnf(ji,jj) + fse3t(ji,jj,jk)
               END DO
            END DO
         END DO
         !
      ELSE IF( ln_rnf_depth_ini ) THEN           ! runoffs applied at the surface
         !
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '    depth of runoff computed once from max value of runoff'
         IF(lwp) WRITE(numout,*) '    max value of the runoff climatologie (over global domain) rn_rnf_max = ', rn_rnf_max
         IF(lwp) WRITE(numout,*) '    depth over which runoffs is spread                        rn_dep_max = ', rn_dep_max
         IF(lwp) WRITE(numout,*) '     create (=1) a runoff depth file or not (=0)      nn_rnf_depth_file  = ', nn_rnf_depth_file

         CALL iom_open( TRIM( sn_rnf%clname ), inum )    !  open runoff file
         CALL iom_gettime( inum, zrec, kntime=nbrec)
         ALLOCATE( zrnfcl(jpi,jpj,nbrec) )     ;      ALLOCATE( zrnf(jpi,jpj) )
         DO jm = 1, nbrec
            CALL iom_get( inum, jpdom_data, TRIM( sn_rnf%clvar ), zrnfcl(:,:,jm), jm )
         END DO
         CALL iom_close( inum )
         zrnf(:,:) = MAXVAL( zrnfcl(:,:,:), DIM=3 )   !  maximum value in time
         DEALLOCATE( zrnfcl )
         !
         h_rnf(:,:) = 1.
         !
         zacoef = rn_dep_max / rn_rnf_max            ! coef of linear relation between runoff and its depth (150m for max of runoff)
         !
         WHERE( zrnf(:,:) > 0._wp )  h_rnf(:,:) = zacoef * zrnf(:,:)   ! compute depth for all runoffs
         !
         DO jj = 1, jpj                     ! take in account min depth of ocean rn_hmin
            DO ji = 1, jpi
               IF( zrnf(ji,jj) > 0._wp ) THEN
                  jk = mbkt(ji,jj)
                  h_rnf(ji,jj) = MIN( h_rnf(ji,jj), gdept_0(ji,jj,jk ) )
               ENDIF
            END DO
         END DO
         !
         nk_rnf(:,:) = 0                       ! number of levels on which runoffs are distributed
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( zrnf(ji,jj) > 0._wp ) THEN
                  jk = 2
                  DO WHILE ( jk /= mbkt(ji,jj) .AND. gdept_0(ji,jj,jk) < h_rnf(ji,jj) ) ;  jk = jk + 1
                  END DO
                  nk_rnf(ji,jj) = jk
               ELSE
                  nk_rnf(ji,jj) = 1
               ENDIF
            END DO
         END DO
         !
         DEALLOCATE( zrnf )
         !
         DO jj = 1, jpj                                ! set the associated depth
            DO ji = 1, jpi
               h_rnf(ji,jj) = 0._wp
               DO jk = 1, nk_rnf(ji,jj)
                  h_rnf(ji,jj) = h_rnf(ji,jj) + fse3t(ji,jj,jk)
               END DO
            END DO
         END DO
         !
         IF( nn_rnf_depth_file == 1 ) THEN      !  save  output nb levels for runoff
            IF(lwp) WRITE(numout,*) '              create runoff depht file'
            CALL iom_open  ( TRIM( sn_dep_rnf%clname ), inum, ldwrt = .TRUE., kiolib = jprstlib )
            CALL iom_rstput( 0, 0, inum, 'rodepth', h_rnf )
            CALL iom_close ( inum )
         ENDIF
      ELSE                                       ! runoffs applied at the surface
         nk_rnf(:,:) = 1
         h_rnf (:,:) = fse3t(:,:,1)
      ENDIF
      !
      rnf(:,:) =  0._wp                         ! runoff initialisation
      rnf_tsc(:,:,:) = 0._wp                    ! runoffs temperature & salinty contents initilisation
      !
      !                                   ! ========================
      !                                   !   River mouth vicinity
      !                                   ! ========================
      !
      IF( ln_rnf_mouth ) THEN                   ! Specific treatment in vicinity of river mouths :
         !                                      !    - Increase Kz in surface layers ( rn_hrnf > 0 )
         !                                      !    - set to zero SSS damping (ln_ssr=T)
         !                                      !    - mixed upstream-centered (ln_traadv_cen2=T)
         !
         IF ( ln_rnf_depth )   CALL ctl_warn( 'sbc_rnf_init: increased mixing turned on but effects may already',   &
            &                                              'be spread through depth by ln_rnf_depth'               )
         !
         nkrnf = 0                                  ! Number of level over which Kz increase
         IF( rn_hrnf > 0._wp ) THEN
            nkrnf = 2
            DO WHILE( nkrnf /= jpkm1 .AND. gdepw_1d(nkrnf+1) < rn_hrnf )   ;   nkrnf = nkrnf + 1
            END DO
            IF( ln_sco )   CALL ctl_warn( 'sbc_rnf: number of levels over which Kz is increased is computed for zco...' )
         ENDIF
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '          Specific treatment used in vicinity of river mouths :'
         IF(lwp) WRITE(numout,*) '             - Increase Kz in surface layers (if rn_hrnf > 0 )'
         IF(lwp) WRITE(numout,*) '               by ', rn_avt_rnf,' m2/s  over ', nkrnf, ' w-levels'
         IF(lwp) WRITE(numout,*) '             - set to zero SSS damping       (if ln_ssr=T)'
         IF(lwp) WRITE(numout,*) '             - mixed upstream-centered       (if ln_traadv_cen2=T)'
         !
         CALL rnf_mouth                             ! set river mouth mask
         !
      ELSE                                      ! No treatment at river mouths
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '          No specific treatment at river mouths'
         rnfmsk  (:,:) = 0._wp
         rnfmsk_z(:)   = 0._wp
         nkrnf = 0
      ENDIF
      !
   END SUBROUTINE sbc_rnf_init


   SUBROUTINE rnf_mouth
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE rnf_mouth  ***
      !!
      !! ** Purpose :   define the river mouths mask
      !!
      !! ** Method  :   read the river mouth mask (=0/1) in the river runoff
      !!                climatological file. Defined a given vertical structure.
      !!                CAUTION, the vertical structure is hard coded on the
      !!                first 5 levels.
      !!                This fields can be used to:
      !!                 - set an upstream advection scheme
      !!                   (ln_rnf_mouth=T and ln_traadv_cen2=T)
      !!                 - increase vertical on the top nn_krnf vertical levels
      !!                   at river runoff input grid point (nn_krnf>=2, see step.F90)
      !!                 - set to zero SSS restoring flux at river mouth grid points
      !!
      !! ** Action  :   rnfmsk   set to 1 at river runoff input, 0 elsewhere
      !!                rnfmsk_z vertical structure
      !!----------------------------------------------------------------------
      INTEGER            ::   inum        ! temporary integers
      CHARACTER(len=255) ::   cl_rnfile   ! runoff file name
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'rnf_mouth : river mouth mask'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~ '
      !
      cl_rnfile = TRIM( cn_dir )//TRIM( sn_cnf%clname )
      IF( .NOT. sn_cnf%ln_clim ) THEN   ;   WRITE(cl_rnfile, '(a,"_y",i4)' ) TRIM( cl_rnfile ), nyear    ! add year
         IF( sn_cnf%cltype == 'monthly' )   WRITE(cl_rnfile, '(a,"m",i2)'  ) TRIM( cl_rnfile ), nmonth   ! add month
      ENDIF
      !
      ! horizontal mask (read in NetCDF file)
      CALL iom_open ( cl_rnfile, inum )                           ! open file
      CALL iom_get  ( inum, jpdom_data, sn_cnf%clvar, rnfmsk )    ! read the river mouth array
      CALL iom_close( inum )                                      ! close file
      !
      IF( nn_closea == 1 )   CALL clo_rnf( rnfmsk )               ! closed sea inflow set as ruver mouth
      !
      rnfmsk_z(:)   = 0._wp                                       ! vertical structure
      rnfmsk_z(1)   = 1.0
      rnfmsk_z(2)   = 1.0                                         ! **********
      rnfmsk_z(3)   = 0.5                                         ! HARD CODED on the 5 first levels
      rnfmsk_z(4)   = 0.25                                        ! **********
      rnfmsk_z(5)   = 0.125
      !
   END SUBROUTINE rnf_mouth

   !!======================================================================
END MODULE sbcrnf
