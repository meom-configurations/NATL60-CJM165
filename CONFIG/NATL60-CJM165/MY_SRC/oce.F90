MODULE oce
   !!======================================================================
   !!                      ***  MODULE  oce  ***
   !! Ocean        :  dynamics and active tracers defined in memory 
   !!======================================================================
   !! History :  1.0  !  2002-11  (G. Madec)  F90: Free form and module
   !!            3.1  !  2009-02  (G. Madec, M. Leclair)  pure z* coordinate
   !!            3.3  !  2010-09  (C. Ethe) TRA-TRC merge: add ts, gtsu, gtsv 4D arrays
   !!----------------------------------------------------------------------
   USE par_oce        ! ocean parameters
   USE lib_mpp        ! MPP library

   IMPLICIT NONE
   PRIVATE

   PUBLIC oce_alloc ! routine called by nemo_init in nemogcm.F90

   LOGICAL, PUBLIC ::   l_traldf_rot = .FALSE.  !: rotated laplacian operator for lateral diffusion

   !! dynamics and tracer fields                            ! before ! now    ! after  ! the after trends becomes the fields
   !! --------------------------                            ! fields ! fields ! trends ! only after tra_zdf and dyn_spg
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   ub   ,  un    , ua     !: i-horizontal velocity          [m/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   vb   ,  vn    , va     !: j-horizontal velocity          [m/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   ua_sv,  va_sv          !: Saved trends (time spliting)   [m/s2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::           wn             !: vertical velocity              [m/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   rotb ,  rotn           !: relative vorticity             [s-1]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   hdivb,  hdivn          !: horizontal divergence          [s-1]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::   tsb  ,  tsn   , tsa    !: 4D T-S fields                  [Celcius,psu] 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::   rab_b,  rab_n          !: thermal/haline expansion coef. [Celcius-1,psu-1]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   rn2b ,  rn2            !: brunt-vaisala frequency**2     [s-2]
   !
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   rhd    !: in situ density anomalie rhd=(rho-rau0)/rau0  [no units]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   rhop   !: potential volumic mass                           [kg/m3]

   !! free surface                                      !  before  ! now    ! after  !
   !! ------------                                      !  fields  ! fields ! fields !
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   ub_b   ,  un_b  ,  ua_b  !: Barotropic velocities at u-point [m/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   vb_b   ,  vn_b  ,  va_b  !: Barotropic velocities at v-point [m/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   sshb   ,  sshn  ,  ssha  !: sea surface height at t-point [m]
   !
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   spgu, spgv               !: horizontal surface pressure gradient

   !! interpolated gradient (only used in zps case)
   !! ---------------------
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   gtsu, gtsv   !: horizontal gradient of T, S bottom u-point 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   gru , grv    !: horizontal gradient of rd at bottom u-point
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   aru , arv    
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   gzu , gzv    
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   ge3ru, ge3rv   !: horizontal gradient of T, S and rd at top v-point  

   !! (ISF) interpolated gradient (only used for ice shelf case) 
   !! --------------------- 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   gtui, gtvi   !: horizontal gradient of T, S and rd at top u-point 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   grui, grvi   !: horizontal gradient of T, S and rd at top v-point  
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   arui, arvi   !: horizontal average  of rd          at top v-point  
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   gzui, gzvi   !: horizontal gradient of z           at top v-point  
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   ge3rui, ge3rvi   !: horizontal gradient of T, S and rd at top v-point  
   !! (ISF) ice load
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   riceload

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   rke          !: kinetic energy

   !! arrays relating to embedding ice in the ocean. These arrays need to be declared 
   !! even if no ice model is required. In the no ice model or traditional levitating 
   !! ice cases they contain only zeros
   !! ---------------------
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   snwice_mass        !: mass of snow and ice at current  ice time step   [Kg/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   snwice_mass_b      !: mass of snow and ice at previous ice time step   [Kg/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   snwice_fmass       !: time evolution of mass of snow+ice               [Kg/m2/s]

   !! Energy budget of the leads (open water embedded in sea ice)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   fraqsr_1lev        !: fraction of solar net radiation absorbed in the first ocean level [-]

   !! Diagnostics for vector Q decomposition, tendancy term of thermal wind imbalance
#if defined key_diavecq   
     REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   ub_vecq   , ua_vecq     !: i-horizontal velocity          [m/s]
     REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   vb_vecq   , va_vecq     !: j-horizontal velocity          [m/s]
#endif
   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO Consortium (2011)
   !! $Id: oce.F90 4990 2014-12-15 16:42:49Z timgraham $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION oce_alloc()
      !!----------------------------------------------------------------------
      !!                   ***  FUNCTION oce_alloc  ***
      !!----------------------------------------------------------------------
      INTEGER :: ierr(5)
      !!----------------------------------------------------------------------
      !
      ALLOCATE( ub   (jpi,jpj,jpk)      , un   (jpi,jpj,jpk)      , ua(jpi,jpj,jpk)       ,     &
         &      vb   (jpi,jpj,jpk)      , vn   (jpi,jpj,jpk)      , va(jpi,jpj,jpk)       ,     &      
         &      ua_sv(jpi,jpj,jpk)      , va_sv(jpi,jpj,jpk)      ,                             &      
         &      wn   (jpi,jpj,jpk)      ,                                                       &
         &      rotb (jpi,jpj,jpk)      , rotn (jpi,jpj,jpk)      ,                             &   
         &      hdivb(jpi,jpj,jpk)      , hdivn(jpi,jpj,jpk)      ,                             &
         &      tsb  (jpi,jpj,jpk,jpts) , tsn  (jpi,jpj,jpk,jpts) , tsa(jpi,jpj,jpk,jpts) ,     &
         &      rab_b(jpi,jpj,jpk,jpts) , rab_n(jpi,jpj,jpk,jpts) ,                             &
         &      rn2b (jpi,jpj,jpk)      , rn2  (jpi,jpj,jpk)                              , STAT=ierr(1) )
         !
      ALLOCATE(rhd (jpi,jpj,jpk) ,                                         &
         &     rhop(jpi,jpj,jpk) ,                                         &
         &     rke(jpi,jpj,jpk)  ,                                         &
         &     sshb(jpi,jpj)     , sshn(jpi,jpj)   , ssha(jpi,jpj)   ,     &
         &     ub_b(jpi,jpj)     , un_b(jpi,jpj)   , ua_b(jpi,jpj)   ,     &
         &     vb_b(jpi,jpj)     , vn_b(jpi,jpj)   , va_b(jpi,jpj)   ,     &
         &     spgu  (jpi,jpj)   , spgv(jpi,jpj)   ,                       &
         &     gtsu(jpi,jpj,jpts), gtsv(jpi,jpj,jpts),                     &
         &     aru(jpi,jpj)      , arv(jpi,jpj)      ,                     &
         &     gzu(jpi,jpj)      , gzv(jpi,jpj)      ,                     &
         &     gru(jpi,jpj)      , grv(jpi,jpj)      ,                     &
         &     ge3ru(jpi,jpj)    , ge3rv(jpi,jpj)    ,                     &
         &     gtui(jpi,jpj,jpts), gtvi(jpi,jpj,jpts),                     &
         &     arui(jpi,jpj)     , arvi(jpi,jpj)     ,                     &
         &     gzui(jpi,jpj)     , gzvi(jpi,jpj)     ,                     &
         &     ge3rui(jpi,jpj)   , ge3rvi(jpi,jpj)   ,                     &
         &     grui(jpi,jpj)     , grvi(jpi,jpj)     ,                     &
         &     riceload(jpi,jpj),                             STAT=ierr(2) )
         !
      ALLOCATE( snwice_mass(jpi,jpj) , snwice_mass_b(jpi,jpj), snwice_fmass(jpi,jpj) , STAT=ierr(3) )
         !
      ALLOCATE( fraqsr_1lev(jpi,jpj) , STAT=ierr(4) )
         !
#if defined key_diavecq   
      ALLOCATE( ub_vecq(jpi,jpj,jpk) , ua_vecq(jpi,jpj,jpk) , &
         &      vb_vecq(jpi,jpj,jpk) , va_vecq(jpi,jpj,jpk) , STAT=ierr(5))
#endif
         !
      oce_alloc = MAXVAL( ierr )
      IF( oce_alloc /= 0 )   CALL ctl_warn('oce_alloc: failed to allocate arrays')
      !
   END FUNCTION oce_alloc

   !!======================================================================
END MODULE oce
