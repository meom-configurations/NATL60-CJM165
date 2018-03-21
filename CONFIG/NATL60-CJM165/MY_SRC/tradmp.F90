MODULE tradmp
   !!======================================================================
   !!                       ***  MODULE  tradmp  ***
   !! Ocean physics: internal restoring trend on active tracers (T and S)
   !!======================================================================
   !! History :  OPA  ! 1991-03  (O. Marti, G. Madec)  Original code
   !!                 ! 1992-06  (M. Imbard)  doctor norme
   !!                 ! 1996-01  (G. Madec)  statement function for e3
   !!                 ! 1997-05  (G. Madec)  macro-tasked on jk-slab
   !!                 ! 1998-07  (M. Imbard, G. Madec) ORCA version
   !!            7.0  ! 2001-02  (M. Imbard)  cofdis, Original code
   !!            8.1  ! 2001-02  (G. Madec, E. Durand)  cleaning
   !!  NEMO      1.0  ! 2002-08  (G. Madec, E. Durand)  free form + modules
   !!            3.2  ! 2009-08  (G. Madec, C. Talandier)  DOCTOR norm for namelist parameter
   !!            3.3  ! 2010-06  (C. Ethe, G. Madec) merge TRA-TRC 
   !!            3.4  ! 2011-04  (G. Madec, C. Ethe) Merge of dtatem and dtasal + suppression of CPP keys
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   tra_dmp_alloc : allocate tradmp arrays
   !!   tra_dmp       : update the tracer trend with the internal damping
   !!   tra_dmp_init  : initialization, namlist read, parameters control
   !!   dtacof        : restoring coefficient for global domain
   !!----------------------------------------------------------------------
   USE oce            ! ocean: variables
   USE dom_oce        ! ocean: domain variables
   USE c1d            ! 1D vertical configuration
   USE trd_oce        ! trends: ocean variables
   USE trdtra         ! trends manager: tracers 
   USE zdf_oce        ! ocean: vertical physics
   USE phycst         ! physical constants
   USE dtatsd         ! data: temperature & salinity
   USE zdfmxl         ! vertical physics: mixed layer depth
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE prtctl         ! Print control
   USE wrk_nemo       ! Memory allocation
   USE timing         ! Timing
   USE iom

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tra_dmp      ! routine called by step.F90
   PUBLIC   tra_dmp_init ! routine called by opa.F90
   PUBLIC   dtacof       ! routine called by tradmp.F90, trcdmp.F90 and dyndmp.F90

   !                                !!* Namelist namtra_dmp : T & S newtonian damping *
   ! nn_zdmp and cn_resto are public as they are used by C1D/dyndmp.F90
   LOGICAL , PUBLIC ::   ln_tradmp   !: internal damping flag
   INTEGER , PUBLIC ::   nn_zdmp     ! = 0/1/2 flag for damping in the mixed layer
   CHARACTER(LEN=200) , PUBLIC :: cn_resto      ! name of netcdf file containing restoration coefficient field
   !

!{ DRAKKAR
   INTEGER , PUBLIC ::   nn_hdmp     ! > 0 standard NEMO CODE
                                     ! = -2 = DRAKKAR customization
                                     ! = -3 = NATL60 customization
   INTEGER , PUBLIC ::   nn_file     ! = 1 create a damping.coeff NetCDF file 
   LOGICAL  ::   ln_dmpmask          !  flag for using mask_dmp file
   INTEGER, SAVE   ::   nk200        !  vertical index for depth > 200 m in ORCA 2
   REAL(wp) ::   rn_timsk            !  restoring time scale used with mask_dmp       [days] 
!}

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   strdmp   !: damping salinity trend (psu/s)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   ttrdmp   !: damping temperature trend (Celcius/s)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   resto    !: restoring coeff. on T and S (s-1)

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: tradmp.F90 4624 2014-04-28 12:09:03Z acc $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION tra_dmp_alloc()
      !!----------------------------------------------------------------------
      !!                ***  FUNCTION tra_dmp_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( strdmp(jpi,jpj,jpk) , ttrdmp(jpi,jpj,jpk), resto(jpi,jpj,jpk), STAT= tra_dmp_alloc )
      !
      IF( lk_mpp            )   CALL mpp_sum ( tra_dmp_alloc )
      IF( tra_dmp_alloc > 0 )   CALL ctl_warn('tra_dmp_alloc: allocation of arrays failed')
      !
   END FUNCTION tra_dmp_alloc


   SUBROUTINE tra_dmp( kt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE tra_dmp  ***
      !!                  
      !! ** Purpose :   Compute the tracer trend due to a newtonian damping
      !!      of the tracer field towards given data field and add it to the
      !!      general tracer trends.
      !!
      !! ** Method  :   Newtonian damping towards t_dta and s_dta computed 
      !!      and add to the general tracer trends:
      !!                     ta = ta + resto * (t_dta - tb)
      !!                     sa = sa + resto * (s_dta - sb)
      !!         The trend is computed either throughout the water column
      !!      (nlmdmp=0) or in area of weak vertical mixing (nlmdmp=1) or
      !!      below the well mixed layer (nlmdmp=2)
      !!
      !! ** Action  : - (ta,sa)   tracer trends updated with the damping trend
      !!----------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      !!
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   zta, zsa             ! local scalars
      REAL(wp), POINTER, DIMENSION(:,:,:,:) ::  zts_dta 
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start( 'tra_dmp')
      !
      CALL wrk_alloc( jpi, jpj, jpk, jpts,  zts_dta )
      !
      !                           !==   input T-S data at kt   ==!
      CALL dta_tsd( kt, zts_dta )            ! read and interpolates T-S data at kt
      !
      SELECT CASE ( nn_zdmp )     !==    type of damping   ==!
      !
      CASE( 0 )                   !==  newtonian damping throughout the water column  ==!
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zta = resto(ji,jj,jk) * ( zts_dta(ji,jj,jk,jp_tem) - tsb(ji,jj,jk,jp_tem) )
                  zsa = resto(ji,jj,jk) * ( zts_dta(ji,jj,jk,jp_sal) - tsb(ji,jj,jk,jp_sal) )
                  tsa(ji,jj,jk,jp_tem) = tsa(ji,jj,jk,jp_tem) + zta
                  tsa(ji,jj,jk,jp_sal) = tsa(ji,jj,jk,jp_sal) + zsa
                  strdmp(ji,jj,jk) = zsa           ! save the trend (used in asmtrj)
                  ttrdmp(ji,jj,jk) = zta      
               END DO
            END DO
         END DO
         !
      CASE ( 1 )                  !==  no damping in the turbocline (avt > 5 cm2/s)  ==!
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  IF( avt(ji,jj,jk) <= 5.e-4_wp ) THEN
                     zta = resto(ji,jj,jk) * ( zts_dta(ji,jj,jk,jp_tem) - tsb(ji,jj,jk,jp_tem) )
                     zsa = resto(ji,jj,jk) * ( zts_dta(ji,jj,jk,jp_sal) - tsb(ji,jj,jk,jp_sal) )
                  ELSE
                     zta = 0._wp
                     zsa = 0._wp  
                  ENDIF
                  tsa(ji,jj,jk,jp_tem) = tsa(ji,jj,jk,jp_tem) + zta
                  tsa(ji,jj,jk,jp_sal) = tsa(ji,jj,jk,jp_sal) + zsa
                  strdmp(ji,jj,jk) = zsa           ! save the salinity trend (used in asmtrj)
                  ttrdmp(ji,jj,jk) = zta
               END DO
            END DO
         END DO
         !
      CASE ( 2 )                  !==  no damping in the mixed layer   ==!
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  IF( fsdept(ji,jj,jk) >= hmlp (ji,jj) ) THEN
                     zta = resto(ji,jj,jk) * ( zts_dta(ji,jj,jk,jp_tem) - tsb(ji,jj,jk,jp_tem) )
                     zsa = resto(ji,jj,jk) * ( zts_dta(ji,jj,jk,jp_sal) - tsb(ji,jj,jk,jp_sal) )
                  ELSE
                     zta = 0._wp
                     zsa = 0._wp  
                  ENDIF
                  tsa(ji,jj,jk,jp_tem) = tsa(ji,jj,jk,jp_tem) + zta
                  tsa(ji,jj,jk,jp_sal) = tsa(ji,jj,jk,jp_sal) + zsa
                  strdmp(ji,jj,jk) = zsa           ! save the salinity trend (used in asmtrj)
                  ttrdmp(ji,jj,jk) = zta
               END DO
            END DO
         END DO
         !
      END SELECT
      !
      IF( l_trdtra )   THEN       ! trend diagnostic
         CALL trd_tra( kt, 'TRA', jp_tem, jptra_dmp, ttrdmp )
         CALL trd_tra( kt, 'TRA', jp_sal, jptra_dmp, strdmp )
      ENDIF
      !                           ! Control print
      IF(ln_ctl)   CALL prt_ctl( tab3d_1=tsa(:,:,:,jp_tem), clinfo1=' dmp  - Ta: ', mask1=tmask,   &
         &                       tab3d_2=tsa(:,:,:,jp_sal), clinfo2=       ' Sa: ', mask2=tmask, clinfo3='tra' )
      !
      CALL wrk_dealloc( jpi, jpj, jpk, jpts,  zts_dta )
      !
      IF( nn_timing == 1 )  CALL timing_stop( 'tra_dmp')
      !
   END SUBROUTINE tra_dmp


   SUBROUTINE tra_dmp_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_dmp_init  ***
      !! 
      !! ** Purpose :   Initialization for the newtonian damping 
      !!
      !! ** Method  :   read the namtra_dmp namelist and check the parameters
      !!----------------------------------------------------------------------
      NAMELIST/namtra_dmp/ ln_tradmp, nn_zdmp, cn_resto
      INTEGER ::  ios         ! Local integer for output status of namelist read
!{DRAKKAR
      INTEGER :: jk ! dummy loop index
      NAMELIST/namtra_dmp/ nn_hdmp , nn_file, ln_dmpmask, rn_timsk
!}
      !!----------------------------------------------------------------------
      !
      REWIND( numnam_ref )   ! Namelist namtra_dmp in reference namelist : T & S relaxation
      READ  ( numnam_ref, namtra_dmp, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namtra_dmp in reference namelist', lwp )
      !
      REWIND( numnam_cfg )   ! Namelist namtra_dmp in configuration namelist : T & S relaxation
      READ  ( numnam_cfg, namtra_dmp, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namtra_dmp in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, namtra_dmp )
      
      IF(lwp) THEN                       ! Namelist print
         WRITE(numout,*)
         WRITE(numout,*) 'tra_dmp_init : T and S newtonian relaxation'
         WRITE(numout,*) '~~~~~~~'
         WRITE(numout,*) '   Namelist namtra_dmp : set relaxation parameters'
         WRITE(numout,*) '      Apply relaxation   or not       ln_tradmp = ', ln_tradmp
         WRITE(numout,*) '      mixed layer damping option      nn_zdmp   = ', nn_zdmp
         WRITE(numout,*) '      Damping file name               cn_resto  = ', cn_resto
!{DRAKKAR
         WRITE(numout,*) '      T and S damping option          nn_hdmp   = ', nn_hdmp
         WRITE(numout,*) '      use a mask_dmp file (T/F)      ln_dmpmask = ', ln_dmpmask
         WRITE(numout,*) '      time scale (mask_dmp)            rn_timsk = ', rn_timsk
         WRITE(numout,*) '      create a damping.coeff file     nn_file   = ', nn_file
!}
         WRITE(numout,*)
      ENDIF

      IF( ln_tradmp) THEN
         !
         !Allocate arrays
         IF( tra_dmp_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'tra_dmp_init: unable to allocate arrays' ) 

         SELECT CASE ( nn_hdmp )
         CASE (  -2  )   ;   IF(lwp) WRITE(numout,*) '   Drakkar customization '
         CASE (  -3  )   ;   IF(lwp) WRITE(numout,*) '   NATL60 customization '
         CASE DEFAULT
                             IF(lwp) WRITE(numout,*) '   Standard Nemo relaxation '
         END SELECT

         !Check values of nn_zdmp
         SELECT CASE ( nn_zdmp )
         CASE ( 0 )  ; IF(lwp) WRITE(numout,*) '   tracer damping as specified by mask'
         CASE ( 1 )  ; IF(lwp) WRITE(numout,*) '   no tracer damping in the turbocline'
         CASE ( 2 )  ; IF(lwp) WRITE(numout,*) '   no tracer damping in the mixed layer'
         END SELECT

         !TG: Initialisation of dtatsd - Would it be better to have dmpdta routine
         !so can damp to something other than intitial conditions files?
         ! JMM : DRAKKAR version is doing so !
         IF( .NOT.ln_tsd_tradmp ) THEN
            CALL ctl_warn( 'tra_dmp_init: read T-S data not initialized, we force ln_tsd_tradmp=T' )
            CALL dta_tsd_init( ld_tradmp=ln_tradmp )        ! forces the initialisation of T-S data
         ENDIF

         !initialise arrays - Are these actually used anywhere else?
         strdmp(:,:,:) = 0._wp
         ttrdmp(:,:,:) = 0._wp

         CALL dtacof( nn_hdmp, nn_file, 'TRA', resto )

         !
         !{ DRAKKAR : for ORCA  damping (used only for ORCA 2 so far ) 
         IF ( cp_cfg == 'orca' ) THEN
           nk200 = jpk
           DO jk = jpk, 1, -1
            IF ( gdept_1d(jk) > 200 ) THEN 
              nk200 = jk
            ELSE
              EXIT
            ENDIF
           ENDDO
         ENDIF
       
      ENDIF
      !
   END SUBROUTINE tra_dmp_init


   SUBROUTINE dtacof( kn_hdmp, kn_file, cdtype , presto )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dtacof  ***
      !!
      !! ** Purpose :   Compute the damping coefficient
      !!
      !! ** Method  :   Arrays defining the damping are computed for each grid
      !!                point for temperature and salinity (resto)
      !!                Damping depends on distance to coast, depth and latitude
      !!
      !! ** Action  : - resto, the damping coeff. for T and S
      !!----------------------------------------------------------------------
      USE ioipsl
      !!
      INTEGER                         , INTENT(in   )  ::  kn_hdmp    ! damping option
      INTEGER                         , INTENT(in   )  ::  kn_file    ! save the damping coef on a file or not
      CHARACTER(len=3)                , INTENT(in   )  ::  cdtype     ! =TRA, TRC or DYN (tracer/dynamics indicator)
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout)  ::  presto     ! restoring coeff. (s-1)
      !
      INTEGER  ::   ji, jj, jk                  ! dummy loop indices
      INTEGER  ::   ii0, ii1, ij0, ij1          ! local integers
      INTEGER  ::   imask        ! File handle 
!{DRAKKAR
      INTEGER  ::   jrelax                      ! width of buffer zone
      INTEGER  ::   inum                        ! Logical unit for reading dmp_mask
      REAL(wp) ::   ztrelax, ztvanish           ! restoring time scale
      REAL(wp) ::   zlon1, zlon2                ! Longitude min and max for patch
      REAL(wp) ::   zbw, zd1, zd2               ! Band width, depth limit
      REAL(wp) ::   zv1, zv2                    ! local scalars
!}
      REAL(wp) ::   zinfl, zlon                 ! local scalars
      REAL(wp) ::   zlat, zlat0, zlat1, zlat2   !   -      -
      REAL(wp) ::   zsdmp, zbdmp                !   -      -
      CHARACTER(len=80)                   :: cfile
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zdct 
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('dtacof')
      !
      CALL wrk_alloc( jpi, jpj, jpk, zdct )

      presto(:,:,:) = 0._wp
      !
      IF( kn_hdmp > 0 ) THEN      ! use standard NEMO code (read from file )
         !initialise arrays - Are these actually used anywhere else?
         strdmp(:,:,:) = 0._wp
         ttrdmp(:,:,:) = 0._wp

         !Read in mask from file
         CALL iom_open ( cn_resto, imask)
         CALL iom_get  ( imask, jpdom_autoglo, 'resto', resto)
         CALL iom_close( imask )

!{ DRAKKAR
      ELSE IF ( nn_hdmp == -2 ) THEN  !
       !--------------------------------------------------------------------------------------------------------------
       ! 3D restoring in semi enclosed seas :
       ! Black Sea :
       ! Lonmin     ! Lonmax       ! Latmin        ! Latmax       ! Band      ! min relax time   !
       !  deg E     ! deg E        !  deg N        ! Deg N        !  deg      !   days           !
       zlon1 = 27.4 ; zlon2 = 42.0 ;  zlat1 = 41.0 ; zlat2 = 47.5 ; zbw = 0.  ; ztrelax = 180. 

       CALL resto_patch ( zlon1, zlon2, zlat1, zlat2, zbw, ztrelax , presto )

       ! Red Sea
       ! Lonmin     ! Lonmax       ! Latmin        ! Latmax       ! Band      ! min relax time   !
       !  deg E     ! deg E        !  deg N        ! Deg N        !  deg      !   days           !
       zlon1 = 29.4 ; zlon2 = 43.6 ;  zlat1 = 12.9 ; zlat2 = 30.3 ; zbw = 0.  ; ztrelax = 180. 

       CALL resto_patch ( zlon1, zlon2, zlat1, zlat2, zbw, ztrelax , presto )

       ! Persian Gulf
       ! Lonmin     ! Lonmax       ! Latmin        ! Latmax       ! Band      ! min relax time   !
       !  deg E     ! deg E        !  deg N        ! Deg N        !  deg      !   days           !
       zlon1 = 46.5 ; zlon2 = 57.0 ;  zlat1 = 23.0 ; zlat2 = 31.5 ; zbw = 1.  ; ztrelax = 180. 

       CALL resto_patch ( zlon1, zlon2, zlat1, zlat2, zbw, ztrelax , presto )
      
       ! ---------------------------------------------------------------------------------------------------------------
       ! 3D overflow restoring 
       ! Gibraltar (limited to 600-1300 m depth range )
       ! Lonmin     ! Lonmax       ! Latmin        ! Latmax       ! Radius    ! relax time   ! dep min    ! dep max    !
       !  deg E     ! deg E        !  deg N        ! Deg N        !  km       !   days       !  m         !   m        !
       zlon1 = -7.  ; zlon2 = -7.  ;  zlat1 = 36.0 ; zlat2 = 36.0 ; zbw = 80. ; ztrelax = 6. ; zd1 = 600. ; zd2 = 1300.

       CALL resto_patch ( zlon1, zlon2, zlat1, zlat2, zbw, ztrelax, presto , zd1, zd2 )

       ! Bab El Mandeb  ( all water column )
       ! Lonmin     ! Lonmax       ! Latmin        ! Latmax       ! Radius    ! relax time   ! dep min    ! dep max    !
       !  deg E     ! deg E        !  deg N        ! Deg N        !  km       !   days       !  m         !   m        !
       zlon1=44.75  ; zlon2 =44.75 ;  zlat1 = 11.5 ; zlat2 = 11.5 ; zbw =100. ; ztrelax = 6. ; zd1 = 0.   ; zd2 = 6000.!

       CALL resto_patch ( zlon1, zlon2, zlat1, zlat2, zbw, ztrelax, presto , zd1, zd2 )

       ! Ormuz Strait  ( all water column )
       ! Lonmin     ! Lonmax       ! Latmin        ! Latmax       ! Radius    ! relax time   ! dep min    ! dep max    !
       !  deg E     ! deg E        !  deg N        ! Deg N        !  km       !   days       !  m         !   m        !
       zlon1=57.75  ; zlon2 =57.75 ;  zlat1 = 25.0 ; zlat2 = 25.0 ; zbw =100. ; ztrelax = 6. ; zd1 = 0.   ; zd2 = 6000.!

       CALL resto_patch ( zlon1, zlon2, zlat1, zlat2, zbw, ztrelax, presto , zd1, zd2 )

       !---------------------------------------------------------------------------------------------------------------
       ! Use a mask_dmp.nc file 
       IF ( ln_dmpmask ) THEN 
         ! in this case, a real 0-1 mask is read into mask_dmp.nc file for 3D restoring following 
         ! particular geometry. This mask is to be build by preprocessing using its own criteria
         ! eg : used for restoring particular water mass. 
         ! The typical restoring time scale is introduced here.
          cfile='dmp_mask.nc'
          CALL iom_open ( cfile, imask )
          CALL iom_get ( imask, jpdom_data, 'wdmp', zdct )  ! use zdct as temporary array
          CALL iom_close (imask)
          WHERE ( zdct > 1 ) zdct = 0.  !  JMM : WHY ???
          ! it seems that this where is not well accepted on brodie => replaced by a loop
          presto(:,:,:) = presto(:,:,:) + zdct(:,:,:)/rn_timsk/86400.
           IF (lwp) WRITE(numout,*) 'dtacof : read dmp_mask.nc file '
           IF (lwp) WRITE(numout,*) '~~~~~'
       ENDIF
       ! Particular cases
       IF ( cp_cfg == 'natl' ) THEN  ! ALL NATL config have an eastern  Med sea truncated
         ! Eastern Med Sea truncated in NATL
         ! Lonmin     ! Lonmax       ! Latmin        ! Latmax       ! Band      ! relax time   !
         !  deg E     ! deg E        !  deg N        ! Deg N        !  deg      !   days       !
         zlon1= 20.0  ; zlon2 = 26.0 ;  zlat1 = 28.0 ; zlat2 = 45.0 ; zbw =2.   ; ztrelax = 3. 

         CALL resto_patch ( zlon1, zlon2, zlat1, zlat2, zbw, ztrelax , presto )
        
!         IF ( Agrif_Root() ) THEN
!         ! Additional modification in case of closed boundary (north and/or south )
!         IF ( .NOT. lk_obc ) THEN
!           jrelax = 28 ! number of point for the buffer zone
!                       ! use resto_patch with lon lat expressed in i,j points, but set ld_ij flag
!           ztrelax  = 3.
!           ztvanish= 100.
!           zv1 = ( ztvanish -  ztrelax ) / ( jrelax  - 1 )
!           DO jj=1, jrelax
!              zv2 = 1./ ( ztrelax + ( jj - 1 ) * zv1 ) / 86400.
!              ! South
!              presto (:, mj0(jj):mj1(jj)                , :) = zv2
!              ! North
!              presto (:, mj0(jpjglo-jj):mj1(jpjglo - jj), :) = zv2
!           ENDDO
!         ENDIF
!         ENDIF  ! agrif root
       ENDIF
      ELSE IF ( nn_hdmp == -3 ) THEN  !   NATL60 configuration
         IF ( cp_cfg == 'natl' .AND. jp_cfg == 60 )  THEN
         ! restoring TS near open boundaries (N/S)
           jrelax = 150 ! number of point for the buffer zone
           ztrelax  = 3.
           ztvanish= 100.
           zv1 = ( ztvanish -  ztrelax ) / ( jrelax  - 1 )
           DO jj=1, jrelax
              zv2 = 1./ ( ztrelax + ( jj - 1 ) * zv1 ) / 86400.
              ! South
              presto (:, mj0(jj):mj1(jj)                , :) = zv2
              ! North
              presto (:, mj0(jpjglo-jj):mj1(jpjglo - jj), :) = zv2
           ENDDO
         ! restoring TS near Eastern boundary
           jrelax = 150 ! number of point for the buffer zone
           ztrelax  = 3.
           ztvanish= 100.
           zv1 = ( ztvanish -  ztrelax ) / ( jrelax  - 1 )
           DO ji=1, jrelax
              zv2 = 1./ ( ztrelax + ( ji - 1 ) * zv1 ) / 86400.
              ! East
!             presto (:, mj0(jpjglo-jj):mj1(jpjglo - jj), :) = zv2
              presto ( mi0(jpiglo-ji):mi1(jpiglo - ji),:, :) = zv2
           ENDDO

         ELSE
            CALL ctl_stop (' You request NATL60 damping for another config ')
         ENDIF

         !                         !--------------------!
      ELSE                         !     No damping     !
         !                         !--------------------!
         CALL ctl_stop( 'Choose a correct value of nn_hdmp or put ln_tradmp to FALSE' )
      ENDIF
      !                            !--------------------------------!
      IF( kn_file == 1 ) THEN      !  save damping coef. in a file  !
         !                         !--------------------------------!
         IF(lwp) WRITE(numout,*) '              create damping.coeff.nc file'
         IF( cdtype == 'TRA' ) cfile = 'damping.coeff'
         IF( cdtype == 'TRC' ) cfile = 'damping.coeff.trc'
         IF( cdtype == 'DYN' ) cfile = 'damping.coeff.dyn'
         cfile = TRIM( cfile )
         CALL iom_open  ( cfile, imask, ldwrt = .TRUE., kiolib = jprstlib )
         CALL iom_rstput( 0, 0, imask, 'Resto', presto )
         CALL iom_close ( imask )
      ENDIF
      !
      CALL wrk_dealloc( jpi, jpj, jpk, zdct )
      !
      IF( nn_timing == 1 )  CALL timing_stop('dtacof')
      !
   END SUBROUTINE dtacof

   SUBROUTINE resto_patch ( plon1, plon2, plat1, plat2, pbw ,ptmax, presto, pz1, pz2 )
      !!------------------------------------------------------------------------
      !!                 ***  Routine resto_patch  ***
      !!
      !! ** Purpose :   modify resto array on a geographically defined zone.
      !!
      !! ** Method  :  Use glamt, gphit arrays. If the defined zone is outside 
      !!              the domain, resto is unchanged. If pz1 and pz2 are provided
      !!              then plon1, plat1 is taken as the position of a the center
      !!              of a circle with decay radius is pbw (in km) 
      !!
      !! ** Action  : IF not present pz1, pz2 : 
      !!              - plon1, plon2 : min and max longitude of the zone (Deg)
      !!              - plat1, plat2 : min and max latitude of the zone (Deg)
      !!              - pbw : band width of the linear decaying restoring (Deg)
      !!              - ptmax : restoring time scale for the inner zone (days)
      !!              IF present pz1 pz2
      !!              - plon1, plat1 : position of the center of the circle
      !!              - pbw = radius (km) of the restoring circle
      !!              - ptmax = time scale at maximum restoring
      !!              - pz1, pz2 : optional: if used, define the depth range (m)
      !!                          for restoring. If not all depths are considered
      !!------------------------------------------------------------------------
      REAL(wp),                   INTENT(in   ) :: plon1, plon2, plat1, plat2, pbw, ptmax
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) :: presto 
      REAL(wp), OPTIONAL, INTENT(in) :: pz1, pz2
      !!
      INTEGER :: ji,jj, jk    ! dummy loop index
      INTEGER :: ik1, ik2     ! limiting vertical index corresponding to pz1,pz2
      INTEGER :: ij0, ij1, iiO, ii1
      INTEGER, DIMENSION(1)           :: iloc 

      REAL(wp) :: zv1, zv2, zv3, zv4, zcoef, ztmp, zdist, zradius2, zcoef2
      REAL(wp), POINTER, DIMENSION(:,:) :: zpatch
      REAL(wp), POINTER, DIMENSION(:)   :: zmask
      !!------------------------------------------------------------------------
      CALL wrk_alloc( jpk,      zmask  )
      CALL wrk_alloc( jpi, jpj, zpatch )
 
      zpatch = 0._wp
      zcoef  = 1._wp/ptmax/86400._wp

      IF (PRESENT (pz1) ) THEN 
        ! horizontal extent
        zradius2 = pbw * pbw !  radius squared
        DO jj = 1, jpj
          DO ji = 1 , jpi
            zpatch(ji,jj) =  sin(gphit(ji,jj)*rad)* sin(plat1*rad)  &
       &                         + cos(gphit(ji,jj)*rad)* cos(plat1*rad)  &
       &                         * cos(rad*(plon1-glamt(ji,jj)))
          ENDDO
        ENDDO

        WHERE ( abs (zpatch ) > 1 ) zpatch = 1.
        DO jj = 1, jpj
          DO ji= 1, jpi 
             ztmp = zpatch(ji,jj)
             zdist = atan(sqrt( (1.-ztmp)/(1+ztmp)) )*2.*ra/1000.
             zpatch(ji,jj) = exp( - zdist*zdist/zradius2 )
          ENDDO
        ENDDO
        ! clean cut off
        WHERE (ABS(zpatch) < 0.01 ) zpatch = 0.
        ! Vertical limitation
        zmask = 1.
        WHERE ( gdept_1d < pz1 .OR. gdept_1d > pz2 ) zmask = 0.
        iloc=MAXLOC(zmask) ; ik1 = iloc(1)
        zmask(1:ik1) = 1.
        iloc=MAXLOC(zmask) ; ik2 = iloc(1) - 1
        zmask(1:ik1) = 1.
        iloc=MAXLOC(zmask) ; ik2 = iloc(1) - 1
        IF (ik2 > 2 ) THEN
          zmask = 0._wp
          zmask(ik1       ) = 0.25_wp
          zmask(ik1+1     ) = 0.75_wp
          zmask(ik1+2:ik2-2) = 1.0_wp
          zmask(ik2-1     ) = 0.75_wp
          zmask(ik2       ) = 0.25_wp
        ELSE
          zmask = 1.   ! all the water column is restored the same
        ENDIF

        DO jk=1, jpk
          presto(:,:,jk)= presto(:,:,jk) + zpatch * zcoef * zmask(jk) 
        ENDDO
        ! JMM : eventually add some checking to avoid locally large resto.

      ELSE
        ! horizontal extent
        zcoef2=1./(pbw +1.e-20 ) ! to avoid division by 0
        DO jj=1,jpj
          DO ji=1,jpi
             zv1=MAX(0., zcoef2*( glamt(ji,jj) - plon1)  )
             zv2=MAX(0., zcoef2*( plon2 - glamt(ji,jj))  )
             zv3=MAX(0., zcoef2*( gphit(ji,jj) - plat1)  )
             zv4=MAX(0., zcoef2*( plat2 - gphit(ji,jj))  )
             zpatch(ji,jj)= MIN( 1., MIN( 1., zv1,zv2,zv3,zv4 ) )
          ENDDO
        ENDDO
          ! resto all over the water column
          presto(:,:,1)= presto(:,:,1) + zpatch *zcoef
          WHERE (zpatch /= 0 ) presto(:,:,1) = MIN( presto(:,:,1), zcoef )
          DO jk=2, jpk
             presto(:,:,jk)=presto(:,:,1)
          ENDDO
      ENDIF
      !
      CALL wrk_dealloc( jpk,       zmask  )
      CALL wrk_dealloc( jpi, jpj , zpatch )
      !
   END SUBROUTINE resto_patch

END MODULE tradmp
