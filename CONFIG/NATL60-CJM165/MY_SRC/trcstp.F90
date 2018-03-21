MODULE trcstp
   !!======================================================================
   !!                       ***  MODULE trcstp  ***
   !! Time-stepping    : time loop of opa for passive tracer
   !!======================================================================
   !! History :  1.0  !  2004-03  (C. Ethe)  Original
   !!----------------------------------------------------------------------
#if defined key_top
   !!----------------------------------------------------------------------
   !!   trc_stp      : passive tracer system time-stepping
   !!----------------------------------------------------------------------
   USE oce_trc          ! ocean dynamics and active tracers variables
   USE sbc_oce
   USE trc
   USE trctrp           ! passive tracers transport
   USE trcsms           ! passive tracers sources and sinks
   USE prtctl_trc       ! Print control for debbuging
   USE trcdia
   USE trcwri
   USE trcrst
   USE trdtrc_oce
   USE trdmxl_trc
   USE iom
   USE in_out_manager
   USE trcsub
   USE diasub


   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_stp    ! called by step

   REAL(wp), DIMENSION(:,:,:), SAVE, ALLOCATABLE ::   qsr_arr ! save qsr during TOP time-step
   REAL(wp) :: rdt_sampl
   INTEGER  :: nb_rec_per_day
   INTEGER  :: isecfst, iseclast
   LOGICAL  :: llnew

   !! * Substitutions
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcstp.F90 6204 2016-01-04 13:47:06Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_stp( kt )
      !!-------------------------------------------------------------------
      !!                     ***  ROUTINE trc_stp  ***
      !!                      
      !! ** Purpose : Time loop of opa for passive tracer
      !! 
      !! ** Method  : 
      !!              Compute the passive tracers trends 
      !!              Update the passive tracers
      !!-------------------------------------------------------------------
      INTEGER, INTENT( in ) ::  kt      ! ocean time-step index
      INTEGER               ::  jk, jn  ! dummy loop indices
      REAL(wp)              ::  ztrai
      CHARACTER (len=25)    ::  charout 

      !!-------------------------------------------------------------------
      !
      IF( nn_timing == 1 )   CALL timing_start('trc_stp')
      !
      IF( kt == nittrc000 .AND. lk_trdmxl_trc )  CALL trd_mxl_trc_init    ! trends: Mixed-layer
      !
      IF( lk_vvl ) THEN                                                   ! update ocean volume due to ssh temporal evolution
         DO jk = 1, jpk
            cvol(:,:,jk) = e1e2t(:,:) * fse3t(:,:,jk) * tmask(:,:,jk)
         END DO
         IF( lk_degrad )  cvol(:,:,:) = cvol(:,:,:) * facvol(:,:,:)       ! degrad option: reduction by facvol
         areatot         = glob_sum( cvol(:,:,:) )
      ENDIF
      !
      IF( l_trcdm2dc )   CALL trc_mean_qsr( kt )
      !   
      IF( lk_diasub )       CALL dia_sub( kt )
      !    
      IF( nn_dttrc /= 1 )   CALL trc_sub_stp( kt )  ! averaging physical variables for sub-stepping
      !    
      IF( MOD( kt , nn_dttrc ) == 0 ) THEN      ! only every nn_dttrc time step
         !
         IF(ln_ctl) THEN
            WRITE(charout,FMT="('kt =', I4,'  d/m/y =',I2,I2,I4)") kt, nday, nmonth, nyear
            CALL prt_ctl_trc_info(charout)
         ENDIF
         !
         tra(:,:,:,:) = 0.e0
         !
                                   CALL trc_rst_opn  ( kt )       ! Open tracer restart file 
         IF( lrst_trc )            CALL trc_rst_cal  ( kt, 'WRITE' )   ! calendar
         IF( lk_iomput ) THEN  ;   CALL trc_wri      ( kt )       ! output of passive tracers with iom I/O manager
         ELSE                  ;   CALL trc_dia      ( kt )       ! output of passive tracers with old I/O manager
         ENDIF
                                   CALL trc_sms      ( kt )       ! tracers: sinks and sources
                                   CALL trc_trp      ( kt )       ! transport of passive tracers
         IF( kt == nittrc000 ) THEN
            CALL iom_close( numrtr )       ! close input tracer restart file
            IF(lwm) CALL FLUSH( numont )   ! flush namelist output
         ENDIF
         IF( lrst_trc )            CALL trc_rst_wri  ( kt )       ! write tracer restart file
         IF( lk_trdmxl_trc  )      CALL trd_mxl_trc  ( kt )       ! trends: Mixed-layer
         !
         IF( nn_dttrc /= 1   )     CALL trc_sub_reset( kt )       ! resetting physical variables when sub-stepping
         !
      ENDIF
      !
!     ztrai = 0._wp                                                   !  content of all tracers
!     DO jn = 1, jptra
!        ztrai = ztrai + glob_sum( trn(:,:,:,jn) * cvol(:,:,:)   )
!     END DO
!     IF( lwp ) WRITE(numstr,9300) kt,  ztrai / areatot
9300  FORMAT(i10,e18.10)
      !
      IF( nn_timing == 1 )   CALL timing_stop('trc_stp')
      !
   END SUBROUTINE trc_stp

   SUBROUTINE trc_mean_qsr( kt )
      !!----------------------------------------------------------------------
      !!             ***  ROUTINE trc_mean_qsr  ***
      !!
      !! ** Purpose :  Compute daily mean qsr for biogeochemical model in case
      !!               of diurnal cycle
      !!
      !! ** Method  : store in TOP the qsr every hour ( or every time-step if the latter 
      !!              is greater than 1 hour ) and then, compute the  mean with 
      !!              a moving average over 24 hours. 
      !!              In coupled mode, the sampling is done at every coupling frequency 
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt
      INTEGER  :: jn

      IF( kt == nittrc000 ) THEN
         IF( ln_cpl )  THEN  
            rdt_sampl = 86400. / ncpl_qsr_freq
            nb_rec_per_day = ncpl_qsr_freq
         ELSE  
            rdt_sampl = MAX( 3600., rdt * nn_dttrc )
            nb_rec_per_day = INT( 86400 / rdt_sampl )
         ENDIF
         !
         IF( lwp ) THEN
            WRITE(numout,*) 
            WRITE(numout,*) ' Sampling frequency dt = ', rdt_sampl, 's','   Number of sampling per day  nrec = ', nb_rec_per_day
            WRITE(numout,*) 
         ENDIF
         !
         !                                            !* Restart: read in restart file
         IF( ln_rsttr .AND. iom_varid( numrtr, 'qsr_mean', ldstop = .FALSE. ) > 0 ) THEN 
            IF(lwp) WRITE(numout,*) 'trc_qsr_mean:   qsr_mean read in the restart file'
            CALL iom_get( numrtr, jpdom_autoglo, 'qsr_mean', qsr_mean )   !  A mean of qsr
         ELSE                                         !* no restart: set from nit000 values
            IF(lwp) WRITE(numout,*) 'trc_qsr_mean:   qsr_mean set to nit000 values'
            qsr_mean(:,:) = qsr(:,:)
         ENDIF
         !
         ALLOCATE( qsr_arr(jpi,jpj,nb_rec_per_day ) )
         DO jn = 1, nb_rec_per_day
             qsr_arr(:,:,jn) = qsr_mean(:,:)
         ENDDO
         !
         isecfst  = nsec_year + nsec1jan000   !   number of seconds between Jan. 1st 00h of nit000 year and the middle of time step
         iseclast = isecfst
         !
      ENDIF
      !
      iseclast = nsec_year + nsec1jan000
      llnew   = ( iseclast - isecfst )  > INT( rdt_sampl )   !   new shortwave to store
      IF( kt /= nittrc000 .AND. llnew ) THEN
          IF( lwp ) WRITE(numout,*) ' New shortwave to sample for TOP at time kt = ', kt, &
             &                      ' time = ', (iseclast+rdt*nn_dttrc/2.)/3600.,'hours '
          isecfst = iseclast
          DO jn = 1, nb_rec_per_day - 1
             qsr_arr(:,:,jn) = qsr_arr(:,:,jn+1)
          ENDDO
          qsr_arr (:,:,nb_rec_per_day) = qsr(:,:)
          qsr_mean(:,:                ) = SUM( qsr_arr(:,:,:), 3 ) / nb_rec_per_day
      ENDIF
      !
      IF( lrst_trc ) THEN    !* Write the mean of qsr in restart file 
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'trc_mean_qsr : write qsr_mean in restart file  kt =', kt
         IF(lwp) WRITE(numout,*) '~~~~~~~'
         CALL iom_rstput( kt, nitrst, numrtw, 'qsr_mean', qsr_mean(:,:) )
      ENDIF
     !
   END SUBROUTINE trc_mean_qsr

#else
   !!----------------------------------------------------------------------
   !!   Default key                                     NO passive tracers
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_stp( kt )        ! Empty routine
      WRITE(*,*) 'trc_stp: You should not have seen this print! error?', kt
   END SUBROUTINE trc_stp
#endif

   !!======================================================================
END MODULE trcstp
