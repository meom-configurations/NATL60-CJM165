MODULE diawri
   !!======================================================================
   !!                     ***  MODULE  diawri  ***
   !! Ocean diagnostics :  write ocean output files
   !!=====================================================================
   !! History :  OPA  ! 1991-03  (M.-A. Foujols)  Original code
   !!            4.0  ! 1991-11  (G. Madec)
   !!                 ! 1992-06  (M. Imbard)  correction restart file
   !!                 ! 1992-07  (M. Imbard)  split into diawri and rstwri
   !!                 ! 1993-03  (M. Imbard)  suppress writibm
   !!                 ! 1998-01  (C. Levy)  NETCDF format using ioipsl INTERFACE
   !!                 ! 1999-02  (E. Guilyardi)  name of netCDF files + variables
   !!            8.2  ! 2000-06  (M. Imbard)  Original code (diabort.F)
   !!   NEMO     1.0  ! 2002-06  (A.Bozec, E. Durand)  Original code (diainit.F)
   !!             -   ! 2002-09  (G. Madec)  F90: Free form and module
   !!             -   ! 2002-12  (G. Madec)  merge of diabort and diainit, F90
   !!                 ! 2005-11  (V. Garnier) Surface pressure gradient organization
   !!            3.2  ! 2008-11  (B. Lemaire) creation from old diawri
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dia_wri       : create the standart output files
   !!   dia_wri_state : create an output NetCDF file for a single instantaeous ocean state and forcing fields
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers 
   USE dom_oce         ! ocean space and time domain
   USE dynadv, ONLY:    ln_dynadv_vec
   USE zdf_oce         ! ocean vertical physics
   USE ldftra_oce      ! ocean active tracers: lateral physics
   USE ldfdyn_oce      ! ocean dynamics: lateral physics
   USE traldf_iso_grif, ONLY : psix_eiv, psiy_eiv
   USE sol_oce         ! solver variables
   USE sbc_oce         ! Surface boundary condition: ocean fields
   USE sbc_ice         ! Surface boundary condition: ice fields
   USE icb_oce         ! Icebergs
   USE icbdia          ! Iceberg budgets
   USE sbcssr          ! restoring term toward SST/SSS climatology
   USE phycst          ! physical constants
   USE zdfmxl          ! mixed layer
   USE dianam          ! build name of file (routine)
   USE zdfddm          ! vertical  physics: double diffusion
   USE diahth          ! thermocline diagnostics
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE in_out_manager  ! I/O manager
   USE diadimg         ! dimg direct access file format output
   USE iom
   USE ioipsl
   USE dynspg_oce, ONLY: un_adv, vn_adv ! barotropic velocities     

#if defined key_lim2
   USE limwri_2 
#elif defined key_lim3
   USE limwri 
#endif
   USE lib_mpp         ! MPP library
   USE timing          ! preformance summary
   USE wrk_nemo        ! working array

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dia_wri                 ! routines called by step.F90
   PUBLIC   dia_wri_state
   PUBLIC   dia_wri_alloc           ! Called by nemogcm module

   INTEGER ::   nid_T, nz_T, nh_T, ndim_T, ndim_hT   ! grid_T file
   INTEGER ::          nb_T              , ndim_bT   ! grid_T file
   INTEGER ::   nid_U, nz_U, nh_U, ndim_U, ndim_hU   ! grid_U file
   INTEGER ::   nid_V, nz_V, nh_V, ndim_V, ndim_hV   ! grid_V file
   INTEGER ::   nid_W, nz_W, nh_W                    ! grid_W file
   INTEGER ::   ndex(1)                              ! ???
   INTEGER, SAVE, ALLOCATABLE, DIMENSION(:) :: ndex_hT, ndex_hU, ndex_hV
   INTEGER, SAVE, ALLOCATABLE, DIMENSION(:) :: ndex_T, ndex_U, ndex_V
   INTEGER, SAVE, ALLOCATABLE, DIMENSION(:) :: ndex_bT

   !! * Substitutions
#  include "zdfddm_substitute.h90"
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION dia_wri_alloc()
      !!----------------------------------------------------------------------
      INTEGER, DIMENSION(2) :: ierr
      !!----------------------------------------------------------------------
      ierr = 0
      ALLOCATE( ndex_hT(jpi*jpj) , ndex_T(jpi*jpj*jpk) ,     &
         &      ndex_hU(jpi*jpj) , ndex_U(jpi*jpj*jpk) ,     &
         &      ndex_hV(jpi*jpj) , ndex_V(jpi*jpj*jpk) , STAT=ierr(1) )
         !
      dia_wri_alloc = MAXVAL(ierr)
      IF( lk_mpp )   CALL mpp_sum( dia_wri_alloc )
      !
  END FUNCTION dia_wri_alloc

#if defined key_dimgout
   !!----------------------------------------------------------------------
   !!   'key_dimgout'                                      DIMG output file
   !!----------------------------------------------------------------------
#   include "diawri_dimg.h90"

#else
   !!----------------------------------------------------------------------
   !!   Default option                                   NetCDF output file
   !!----------------------------------------------------------------------
# if defined key_iomput
   !!----------------------------------------------------------------------
   !!   'key_iomput'                                        use IOM library
   !!----------------------------------------------------------------------

   SUBROUTINE dia_wri( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dia_wri  ***
      !!                   
      !! ** Purpose :   Standard output of opa: dynamics and tracer fields 
      !!      NETCDF format is used by default 
      !!
      !! ** Method  :  use iom_put
      !!----------------------------------------------------------------------
      !!
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
      !!
      INTEGER                      ::   ji, jj, jk              ! dummy loop indices
      INTEGER                      ::   jkbot                   !
      REAL(wp)                     ::   zztmp, zztmpx, zztmpy   ! 
      !!
      REAL(wp), POINTER, DIMENSION(:,:)   :: z2d       ! 2D workspace
      REAL(wp), POINTER, DIMENSION(:,:,:) :: z3d      ! 3D workspace
      !!----------------------------------------------------------------------
      ! 
      IF( nn_timing == 1 )   CALL timing_start('dia_wri')
      ! 
      CALL wrk_alloc( jpi , jpj      , z2d )
      CALL wrk_alloc( jpi , jpj, jpk , z3d )
      !
      ! Output the initial state and forcings
      IF( ninist == 1 ) THEN                       
         CALL dia_wri_state( 'output.init', kt )
         ninist = 0
      ENDIF

      ! Output of initial vertical scale factor
      CALL iom_put("e3t_0", e3t_0(:,:,:) )
      CALL iom_put("e3u_0", e3t_0(:,:,:) )
      CALL iom_put("e3v_0", e3t_0(:,:,:) )
      !
         CALL iom_put( "e3t" , fse3t_n(:,:,:) )
         CALL iom_put( "e3u" , fse3u_n(:,:,:) )
         CALL iom_put( "e3v" , fse3v_n(:,:,:) )
         CALL iom_put( "e3w" , fse3w_n(:,:,:) )
      IF( iom_use("e3tdef") )   &
         CALL iom_put( "e3tdef"  , ( ( fse3t_n(:,:,:) - e3t_0(:,:,:) ) / e3t_0(:,:,:) * 100 * tmask(:,:,:) ) ** 2 )


      CALL iom_put( "ssh" , sshn )                 ! sea surface height
      
      CALL iom_put( "toce", tsn(:,:,:,jp_tem) )    ! 3D temperature
      CALL iom_put(  "sst", tsn(:,:,1,jp_tem) )    ! surface temperature
      IF ( iom_use("sbt") ) THEN
         DO jj = 1, jpj
            DO ji = 1, jpi
               jkbot = mbkt(ji,jj)
               z2d(ji,jj) = tsn(ji,jj,jkbot,jp_tem)
            END DO
         END DO  
         CALL iom_put( "sbt", z2d )                ! bottom temperature
      ENDIF
      
      CALL iom_put( "soce", tsn(:,:,:,jp_sal) )    ! 3D salinity
      CALL iom_put(  "sss", tsn(:,:,1,jp_sal) )    ! surface salinity
      IF ( iom_use("sbs") ) THEN
         DO jj = 1, jpj
            DO ji = 1, jpi
               jkbot = mbkt(ji,jj)
               z2d(ji,jj) = tsn(ji,jj,jkbot,jp_sal)
            END DO
         END DO  
         CALL iom_put( "sbs", z2d )                ! bottom salinity
      ENDIF

      IF ( iom_use("taubot") ) THEN                ! bottom stress
         z2d(:,:) = 0._wp
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zztmpx = (  bfrua(ji  ,jj) * un(ji  ,jj,mbku(ji  ,jj))  &
                      &  + bfrua(ji-1,jj) * un(ji-1,jj,mbku(ji-1,jj))  )      
               zztmpy = (  bfrva(ji,  jj) * vn(ji,jj  ,mbkv(ji,jj  ))  &
                      &  + bfrva(ji,jj-1) * vn(ji,jj-1,mbkv(ji,jj-1))  ) 
               z2d(ji,jj) = rau0 * SQRT( zztmpx * zztmpx + zztmpy * zztmpy ) * tmask(ji,jj,1) 
               !
            END DO
         END DO  
         CALL lbc_lnk( z2d, 'T', 1. )
         CALL iom_put( "taubot", z2d )           
         ENDIF
         
      CALL iom_put( "uoce", un(:,:,:)         )    ! 3D i-current
      CALL iom_put(  "ssu", un(:,:,1)         )    ! surface i-current
! NATL60
! JMM add subsurface daily output
      CALL iom_put( "uoce10"  , un(:,:,10)                             )    ! surface i-current      
      CALL iom_put( "uoce30"  , un(:,:,15)                             )    ! surface i-current      
      CALL iom_put( "uoce50"  , un(:,:,21)                             )    ! surface i-current      
!
      IF ( iom_use("sbu") ) THEN
            DO jj = 1, jpj
               DO ji = 1, jpi
               jkbot = mbku(ji,jj)
               z2d(ji,jj) = un(ji,jj,jkbot)
               END DO
            END DO
         CALL iom_put( "sbu", z2d )                ! bottom i-current
         ENDIF
#if defined key_dynspg_ts
      CALL iom_put(  "ubar", un_adv(:,:)      )    ! barotropic i-current
#else
      CALL iom_put(  "ubar", un_b(:,:)        )    ! barotropic i-current
#endif
      
      CALL iom_put( "voce", vn(:,:,:)         )    ! 3D j-current
      CALL iom_put(  "ssv", vn(:,:,1)         )    ! surface j-current
! NATL60 
! JMM add subsurface daily output
      CALL iom_put( "voce10"  , vn(:,:,10)                             )    ! surface i-current      
      CALL iom_put( "voce30"  , vn(:,:,15)                             )    ! surface i-current      
      CALL iom_put( "voce50"  , vn(:,:,21)                             )    ! surface i-current      
!
      IF ( iom_use("sbv") ) THEN
            DO jj = 1, jpj
               DO ji = 1, jpi
               jkbot = mbkv(ji,jj)
               z2d(ji,jj) = vn(ji,jj,jkbot)
               END DO
            END DO
         CALL iom_put( "sbv", z2d )                ! bottom j-current
         ENDIF
#if defined key_dynspg_ts
      CALL iom_put(  "vbar", vn_adv(:,:)      )    ! barotropic j-current
#else
      CALL iom_put(  "vbar", vn_b(:,:)        )    ! barotropic j-current
#endif

      CALL iom_put( "woce", wn )                   ! vertical velocity
      IF( iom_use('w_masstr') .OR. iom_use('w_masstr2') ) THEN   ! vertical mass transport & its square value
         ! Caution: in the VVL case, it only correponds to the baroclinic mass transport.
         z2d(:,:) = rau0 * e12t(:,:)
         DO jk = 1, jpk
            z3d(:,:,jk) = wn(:,:,jk) * z2d(:,:)
               END DO
         CALL iom_put( "w_masstr" , z3d )  
         IF( iom_use('w_masstr2') )   CALL iom_put( "w_masstr2", z3d(:,:,:) * z3d(:,:,:) )
         ENDIF

      CALL iom_put( "avt"    , avt                                   )    ! T vert. eddy diff. coef.
      CALL iom_put( "avmu"    , avmu                                  )    ! U vert. eddy visc. coef.
      CALL iom_put( "avmv"    , avmv                                  )    ! U vert. eddy visc. coef.
      CALL iom_put( "avs" , fsavs(:,:,:)               )    ! S vert. eddy diff. coef. (useful only with key_zdfddm)
                                                            ! Log of eddy diff coef
      IF( iom_use('logavt') )   CALL iom_put( "logavt", LOG( MAX( 1.e-20_wp, avt  (:,:,:) ) ) )
      IF( iom_use('logavs') )   CALL iom_put( "logavs", LOG( MAX( 1.e-20_wp, fsavs(:,:,:) ) ) )

      IF ( iom_use("sstgrad") .OR. iom_use("sstgrad2") ) THEN
      DO jj = 2, jpjm1                                    ! sst gradient
         DO ji = fs_2, fs_jpim1   ! vector opt.
            zztmp      = tsn(ji,jj,1,jp_tem)
            zztmpx     = ( tsn(ji+1,jj  ,1,jp_tem) - zztmp ) / e1u(ji,jj) + ( zztmp - tsn(ji-1,jj  ,1,jp_tem) ) / e1u(ji-1,jj  )
            zztmpy     = ( tsn(ji  ,jj+1,1,jp_tem) - zztmp ) / e2v(ji,jj) + ( zztmp - tsn(ji  ,jj-1,1,jp_tem) ) / e2v(ji  ,jj-1)
            z2d(ji,jj) = 0.25 * ( zztmpx * zztmpx + zztmpy * zztmpy )   &
               &              * umask(ji,jj,1) * umask(ji-1,jj,1) * vmask(ji,jj,1) * umask(ji,jj-1,1)
         END DO
      END DO
      CALL lbc_lnk( z2d, 'T', 1. )
      CALL iom_put( "sstgrad2",  z2d               )    ! square of module of sst gradient
      z2d(:,:) = SQRT( z2d(:,:) )
      CALL iom_put( "sstgrad" ,  z2d               )    ! module of sst gradient
      ENDIF

      ! clem: heat and salt content
      IF( iom_use("heatc") ) THEN
      z2d(:,:)  = 0._wp 
      DO jk = 1, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  z2d(ji,jj) = z2d(ji,jj) + fse3t(ji,jj,jk) * tsn(ji,jj,jk,jp_tem) * tmask(ji,jj,jk)
            END DO
         END DO
      END DO
         CALL iom_put( "heatc", (rau0 * rcp) * z2d )    ! vertically integrated heat content (J/m2)
      ENDIF

      IF( iom_use("saltc") ) THEN
         z2d(:,:)  = 0._wp 
         DO jk = 1, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  z2d(ji,jj) = z2d(ji,jj) + fse3t(ji,jj,jk) * tsn(ji,jj,jk,jp_sal) * tmask(ji,jj,jk)
               END DO
            END DO
         END DO
         CALL iom_put( "saltc", rau0 * z2d )   ! vertically integrated salt content (PSU*kg/m2)
      ENDIF
      !
      IF ( iom_use("eken") ) THEN
         rke(:,:,jk) = 0._wp                               !      kinetic energy 
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zztmp   = 1._wp / ( e1e2t(ji,jj) * fse3t(ji,jj,jk) )
                  zztmpx  = 0.5 * (  un(ji-1,jj,jk) * un(ji-1,jj,jk) * e2u(ji-1,jj) * fse3u(ji-1,jj,jk)    &
                     &             + un(ji  ,jj,jk) * un(ji  ,jj,jk) * e2u(ji  ,jj) * fse3u(ji  ,jj,jk) )  &
                     &          *  zztmp 
                  !
                  zztmpy  = 0.5 * (  vn(ji,jj-1,jk) * vn(ji,jj-1,jk) * e1v(ji,jj-1) * fse3v(ji,jj-1,jk)    &
                     &             + vn(ji,jj  ,jk) * vn(ji,jj  ,jk) * e1v(ji,jj  ) * fse3v(ji,jj  ,jk) )  &
                     &          *  zztmp 
                  !
                  rke(ji,jj,jk) = 0.5_wp * ( zztmpx + zztmpy )
                  !
               ENDDO
            ENDDO
         ENDDO
         CALL lbc_lnk( rke, 'T', 1. )
         CALL iom_put( "eken", rke )           
      ENDIF
      !
      CALL iom_put( "hdiv", hdivn )                  ! Horizontal divergence
      !
      IF( iom_use("u_masstr") .OR. iom_use("u_heattr") .OR. iom_use("u_salttr") ) THEN
         z3d(:,:,jpk) = 0.e0
         DO jk = 1, jpkm1
            z3d(:,:,jk) = rau0 * un(:,:,jk) * e2u(:,:) * fse3u(:,:,jk) * umask(:,:,jk)
         END DO
         CALL iom_put( "u_masstr", z3d )                  ! mass transport in i-direction
      ENDIF

      IF( iom_use("u_heattr") ) THEN
         z2d(:,:) = 0.e0 
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  z2d(ji,jj) = z2d(ji,jj) + z3d(ji,jj,jk) * ( tsn(ji,jj,jk,jp_tem) + tsn(ji+1,jj,jk,jp_tem) )
               END DO
            END DO
         END DO
         CALL lbc_lnk( z2d, 'U', -1. )
         CALL iom_put( "u_heattr", (0.5 * rcp) * z2d )    ! heat transport in i-direction
      ENDIF

      IF( iom_use("u_salttr") ) THEN
         z2d(:,:) = 0.e0 
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  z2d(ji,jj) = z2d(ji,jj) + z3d(ji,jj,jk) * ( tsn(ji,jj,jk,jp_sal) + tsn(ji+1,jj,jk,jp_sal) )
               END DO
            END DO
         END DO
         CALL lbc_lnk( z2d, 'U', -1. )
         CALL iom_put( "u_salttr", 0.5 * z2d )            ! heat transport in i-direction
      ENDIF

      
      IF( iom_use("v_masstr") .OR. iom_use("v_heattr") .OR. iom_use("v_salttr") ) THEN
         z3d(:,:,jpk) = 0.e0
         DO jk = 1, jpkm1
            z3d(:,:,jk) = rau0 * vn(:,:,jk) * e1v(:,:) * fse3v(:,:,jk) * vmask(:,:,jk)
         END DO
         CALL iom_put( "v_masstr", z3d )                  ! mass transport in j-direction
      ENDIF

      IF( iom_use("v_heattr") ) THEN
         z2d(:,:) = 0.e0 
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  z2d(ji,jj) = z2d(ji,jj) + z3d(ji,jj,jk) * ( tsn(ji,jj,jk,jp_tem) + tsn(ji,jj+1,jk,jp_tem) )
               END DO
            END DO
         END DO
         CALL lbc_lnk( z2d, 'V', -1. )
         CALL iom_put( "v_heattr", (0.5 * rcp) * z2d )    !  heat transport in j-direction
      ENDIF

      IF( iom_use("v_salttr") ) THEN
         z2d(:,:) = 0.e0 
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  z2d(ji,jj) = z2d(ji,jj) + z3d(ji,jj,jk) * ( tsn(ji,jj,jk,jp_sal) + tsn(ji,jj+1,jk,jp_sal) )
               END DO
            END DO
         END DO
         CALL lbc_lnk( z2d, 'V', -1. )
         CALL iom_put( "v_salttr", 0.5 * z2d )            !  heat transport in j-direction
      ENDIF
      !
      CALL wrk_dealloc( jpi , jpj      , z2d )
      CALL wrk_dealloc( jpi , jpj, jpk , z3d )
      !
      IF( nn_timing == 1 )   CALL timing_stop('dia_wri')
      !
   END SUBROUTINE dia_wri

#else
   !!----------------------------------------------------------------------
   !!   Default option                                  use IOIPSL  library
   !!----------------------------------------------------------------------

   SUBROUTINE dia_wri( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dia_wri  ***
      !!                   
      !! ** Purpose :   Standard output of opa: dynamics and tracer fields 
      !!      NETCDF format is used by default 
      !!
      !! ** Method  :   At the beginning of the first time step (nit000), 
      !!      define all the NETCDF files and fields
      !!      At each time step call histdef to compute the mean if ncessary
      !!      Each nwrite time step, output the instantaneous or mean fields
      !!----------------------------------------------------------------------
      !!
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
      !!
      LOGICAL ::   ll_print = .FALSE.                        ! =T print and flush numout
      CHARACTER (len=80) ::   clhstnam, clop, clmx           ! local names
      INTEGER  ::   inum = 11                                ! temporary logical unit
      INTEGER  ::   ji, jj, jk                               ! dummy loop indices
      INTEGER  ::   ierr                                     ! error code return from allocation
      INTEGER  ::   iimi, iima, ipk, it, itmod, ijmi, ijma   ! local integers
      INTEGER  ::   jn, ierror                               ! local integers
      REAL(wp) ::   zsto, zout, zmax, zjulian, zdt           ! local scalars
      !!
      REAL(wp), POINTER, DIMENSION(:,:)   :: zw2d       ! 2D workspace
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zw3d       ! 3D workspace
      !!----------------------------------------------------------------------
      ! 
      IF( nn_timing == 1 )   CALL timing_start('dia_wri')
      !
      CALL wrk_alloc( jpi , jpj      , zw2d )
      IF ( ln_traldf_gdia .OR. lk_vvl )  call wrk_alloc( jpi , jpj , jpk  , zw3d )
      !
      ! Output the initial state and forcings
      IF( ninist == 1 ) THEN                       
         CALL dia_wri_state( 'output.init', kt )
         ninist = 0
      ENDIF
      !
      ! 0. Initialisation
      ! -----------------

      ! local variable for debugging
      ll_print = .FALSE.
      ll_print = ll_print .AND. lwp

      ! Define frequency of output and means
      zdt = rdt
      IF( nacc == 1 ) zdt = rdtmin
      clop = "x"         ! no use of the mask value (require less cpu time, and otherwise the model crashes)
#if defined key_diainstant
      zsto = nwrite * zdt
      clop = "inst("//TRIM(clop)//")"
#else
      zsto=zdt
      clop = "ave("//TRIM(clop)//")"
#endif
      zout = nwrite * zdt
      zmax = ( nitend - nit000 + 1 ) * zdt

      ! Define indices of the horizontal output zoom and vertical limit storage
      iimi = 1      ;      iima = jpi
      ijmi = 1      ;      ijma = jpj
      ipk = jpk

      ! define time axis
      it = kt
      itmod = kt - nit000 + 1


      ! 1. Define NETCDF files and fields at beginning of first time step
      ! -----------------------------------------------------------------

      IF( kt == nit000 ) THEN

         ! Define the NETCDF files (one per grid)

         ! Compute julian date from starting date of the run
         CALL ymds2ju( nyear, nmonth, nday, rdt, zjulian )
         zjulian = zjulian - adatrj   !   set calendar origin to the beginning of the experiment
         IF(lwp)WRITE(numout,*)
         IF(lwp)WRITE(numout,*) 'Date 0 used :', nit000, ' YEAR ', nyear,   &
            &                    ' MONTH ', nmonth, ' DAY ', nday, 'Julian day : ', zjulian
         IF(lwp)WRITE(numout,*) ' indexes of zoom = ', iimi, iima, ijmi, ijma,   &
                                 ' limit storage in depth = ', ipk

         ! WRITE root name in date.file for use by postpro
         IF(lwp) THEN
            CALL dia_nam( clhstnam, nwrite,' ' )
            CALL ctl_opn( inum, 'date.file', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, numout, lwp, narea )
            WRITE(inum,*) clhstnam
            CLOSE(inum)
         ENDIF

         ! Define the T grid FILE ( nid_T )

         CALL dia_nam( clhstnam, nwrite, 'grid_T' )
         IF(lwp) WRITE(numout,*) " Name of NETCDF file ", clhstnam    ! filename
         CALL histbeg( clhstnam, jpi, glamt, jpj, gphit,           &  ! Horizontal grid: glamt and gphit
            &          iimi, iima-iimi+1, ijmi, ijma-ijmi+1,       &
            &          nit000-1, zjulian, zdt, nh_T, nid_T, domain_id=nidom, snc4chunks=snc4set )
         CALL histvert( nid_T, "deptht", "Vertical T levels",      &  ! Vertical grid: gdept
            &           "m", ipk, gdept_1d, nz_T, "down" )
         !                                                            ! Index of ocean points
         CALL wheneq( jpi*jpj*ipk, tmask, 1, 1., ndex_T , ndim_T  )      ! volume
         CALL wheneq( jpi*jpj    , tmask, 1, 1., ndex_hT, ndim_hT )      ! surface
         !
         IF( ln_icebergs ) THEN
            !
            !! allocation cant go in dia_wri_alloc because ln_icebergs is only set after 
            !! that routine is called from nemogcm, so do it here immediately before its needed
            ALLOCATE( ndex_bT(jpi*jpj*nclasses), STAT=ierror )
            IF( lk_mpp )   CALL mpp_sum( ierror )
            IF( ierror /= 0 ) THEN
               CALL ctl_stop('dia_wri: failed to allocate iceberg diagnostic array')
               RETURN
            ENDIF
            !
            !! iceberg vertical coordinate is class number
            CALL histvert( nid_T, "class", "Iceberg class",      &  ! Vertical grid: class
               &           "number", nclasses, class_num, nb_T )
            !
            !! each class just needs the surface index pattern
            ndim_bT = 3
            DO jn = 1,nclasses
               ndex_bT((jn-1)*jpi*jpj+1:jn*jpi*jpj) = ndex_hT(1:jpi*jpj)
            ENDDO
            !
         ENDIF

         ! Define the U grid FILE ( nid_U )

         CALL dia_nam( clhstnam, nwrite, 'grid_U' )
         IF(lwp) WRITE(numout,*) " Name of NETCDF file ", clhstnam    ! filename
         CALL histbeg( clhstnam, jpi, glamu, jpj, gphiu,           &  ! Horizontal grid: glamu and gphiu
            &          iimi, iima-iimi+1, ijmi, ijma-ijmi+1,       &
            &          nit000-1, zjulian, zdt, nh_U, nid_U, domain_id=nidom, snc4chunks=snc4set )
         CALL histvert( nid_U, "depthu", "Vertical U levels",      &  ! Vertical grid: gdept
            &           "m", ipk, gdept_1d, nz_U, "down" )
         !                                                            ! Index of ocean points
         CALL wheneq( jpi*jpj*ipk, umask, 1, 1., ndex_U , ndim_U  )      ! volume
         CALL wheneq( jpi*jpj    , umask, 1, 1., ndex_hU, ndim_hU )      ! surface

         ! Define the V grid FILE ( nid_V )

         CALL dia_nam( clhstnam, nwrite, 'grid_V' )                   ! filename
         IF(lwp) WRITE(numout,*) " Name of NETCDF file ", clhstnam
         CALL histbeg( clhstnam, jpi, glamv, jpj, gphiv,           &  ! Horizontal grid: glamv and gphiv
            &          iimi, iima-iimi+1, ijmi, ijma-ijmi+1,       &
            &          nit000-1, zjulian, zdt, nh_V, nid_V, domain_id=nidom, snc4chunks=snc4set )
         CALL histvert( nid_V, "depthv", "Vertical V levels",      &  ! Vertical grid : gdept
            &          "m", ipk, gdept_1d, nz_V, "down" )
         !                                                            ! Index of ocean points
         CALL wheneq( jpi*jpj*ipk, vmask, 1, 1., ndex_V , ndim_V  )      ! volume
         CALL wheneq( jpi*jpj    , vmask, 1, 1., ndex_hV, ndim_hV )      ! surface

         ! Define the W grid FILE ( nid_W )

         CALL dia_nam( clhstnam, nwrite, 'grid_W' )                   ! filename
         IF(lwp) WRITE(numout,*) " Name of NETCDF file ", clhstnam
         CALL histbeg( clhstnam, jpi, glamt, jpj, gphit,           &  ! Horizontal grid: glamt and gphit
            &          iimi, iima-iimi+1, ijmi, ijma-ijmi+1,       &
            &          nit000-1, zjulian, zdt, nh_W, nid_W, domain_id=nidom, snc4chunks=snc4set )
         CALL histvert( nid_W, "depthw", "Vertical W levels",      &  ! Vertical grid: gdepw
            &          "m", ipk, gdepw_1d, nz_W, "down" )


         ! Declare all the output fields as NETCDF variables

         !                                                                                      !!! nid_T : 3D
         CALL histdef( nid_T, "votemper", "Temperature"                        , "C"      ,   &  ! tn
            &          jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout )
         CALL histdef( nid_T, "vosaline", "Salinity"                           , "PSU"    ,   &  ! sn
            &          jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout )
         IF(  lk_vvl  ) THEN
            CALL histdef( nid_T, "vovvle3t", "Level thickness"                    , "m"      ,&  ! e3t_n
            &             jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout )
            CALL histdef( nid_T, "vovvldep", "T point depth"                      , "m"      ,&  ! e3t_n
            &             jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout )
            CALL histdef( nid_T, "vovvldef", "Squared level deformation"          , "%^2"    ,&  ! e3t_n
            &             jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout )
         ENDIF
         !                                                                                      !!! nid_T : 2D
         CALL histdef( nid_T, "sosstsst", "Sea Surface temperature"            , "C"      ,   &  ! sst
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         CALL histdef( nid_T, "sosaline", "Sea Surface Salinity"               , "PSU"    ,   &  ! sss
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         CALL histdef( nid_T, "sossheig", "Sea Surface Height"                 , "m"      ,   &  ! ssh
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         CALL histdef( nid_T, "sowaflup", "Net Upward Water Flux"              , "Kg/m2/s",   &  ! (emp-rnf)
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         CALL histdef( nid_T, "sorunoff", "River runoffs"                      , "Kg/m2/s",   &  ! runoffs
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         CALL histdef( nid_T, "sosfldow", "downward salt flux"                 , "PSU/m2/s",  &  ! sfx
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         IF(  .NOT. lk_vvl  ) THEN
         CALL histdef( nid_T, "sosst_cd", "Concentration/Dilution term on temperature"        &  ! emp * tsn(:,:,1,jp_tem)
            &                                                                  , "KgC/m2/s",  &  ! sosst_cd
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         CALL histdef( nid_T, "sosss_cd", "Concentration/Dilution term on salinity"           &  ! emp * tsn(:,:,1,jp_sal)
            &                                                                  , "KgPSU/m2/s",&  ! sosss_cd
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         ENDIF
         CALL histdef( nid_T, "sohefldo", "Net Downward Heat Flux"             , "W/m2"   ,   &  ! qns + qsr
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         CALL histdef( nid_T, "soshfldo", "Shortwave Radiation"                , "W/m2"   ,   &  ! qsr
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         CALL histdef( nid_T, "somixhgt", "Turbocline Depth"                   , "m"      ,   &  ! hmld
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         CALL histdef( nid_T, "somxl010", "Mixed Layer Depth 0.01"             , "m"      ,   &  ! hmlp
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         CALL histdef( nid_T, "soicecov", "Ice fraction"                       , "[0,1]"  ,   &  ! fr_i
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         CALL histdef( nid_T, "sowindsp", "wind speed at 10m"                  , "m/s"    ,   &  ! wndm
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
!
         IF( ln_icebergs ) THEN
            CALL histdef( nid_T, "calving"             , "calving mass input"                       , "kg/s"   , &
               &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
            CALL histdef( nid_T, "calving_heat"        , "calving heat flux"                        , "XXXX"   , &
               &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
            CALL histdef( nid_T, "berg_floating_melt"  , "Melt rate of icebergs + bits"             , "kg/m2/s", &
               &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
            CALL histdef( nid_T, "berg_stored_ice"     , "Accumulated ice mass by class"            , "kg"     , &
               &          jpi, jpj, nh_T, nclasses  , 1, nclasses  , nb_T , 32, clop, zsto, zout )
            IF( ln_bergdia ) THEN
               CALL histdef( nid_T, "berg_melt"           , "Melt rate of icebergs"                    , "kg/m2/s", &
                  &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
               CALL histdef( nid_T, "berg_buoy_melt"      , "Buoyancy component of iceberg melt rate"  , "kg/m2/s", &
                  &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
               CALL histdef( nid_T, "berg_eros_melt"      , "Erosion component of iceberg melt rate"   , "kg/m2/s", &
                  &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
               CALL histdef( nid_T, "berg_conv_melt"      , "Convective component of iceberg melt rate", "kg/m2/s", &
                  &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
               CALL histdef( nid_T, "berg_virtual_area"   , "Virtual coverage by icebergs"             , "m2"     , &
                  &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
               CALL histdef( nid_T, "bits_src"           , "Mass source of bergy bits"                , "kg/m2/s", &
                  &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
               CALL histdef( nid_T, "bits_melt"          , "Melt rate of bergy bits"                  , "kg/m2/s", &
                  &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
               CALL histdef( nid_T, "bits_mass"          , "Bergy bit density field"                  , "kg/m2"  , &
                  &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
               CALL histdef( nid_T, "berg_mass"           , "Iceberg density field"                    , "kg/m2"  , &
                  &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
               CALL histdef( nid_T, "berg_real_calving"   , "Calving into iceberg class"               , "kg/s"   , &
                  &          jpi, jpj, nh_T, nclasses  , 1, nclasses  , nb_T , 32, clop, zsto, zout )
            ENDIF
         ENDIF

         IF( .NOT. ln_cpl ) THEN
         CALL histdef( nid_T, "sohefldp", "Surface Heat Flux: Damping"         , "W/m2"   ,   &  ! qrp
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         CALL histdef( nid_T, "sowafldp", "Surface Water Flux: Damping"        , "Kg/m2/s",   &  ! erp
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         CALL histdef( nid_T, "sosafldp", "Surface salt flux: damping"         , "Kg/m2/s",   &  ! erp * sn
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         ENDIF

         IF( ln_cpl .AND. nn_ice <= 1 ) THEN
         CALL histdef( nid_T, "sohefldp", "Surface Heat Flux: Damping"         , "W/m2"   ,   &  ! qrp
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         CALL histdef( nid_T, "sowafldp", "Surface Water Flux: Damping"        , "Kg/m2/s",   &  ! erp
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         CALL histdef( nid_T, "sosafldp", "Surface salt flux: Damping"         , "Kg/m2/s",   &  ! erp * sn
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         ENDIF
         
         clmx ="l_max(only(x))"    ! max index on a period
         CALL histdef( nid_T, "sobowlin", "Bowl Index"                         , "W-point",   &  ! bowl INDEX 
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clmx, zsto, zout )
#if defined key_diahth
         CALL histdef( nid_T, "sothedep", "Thermocline Depth"                  , "m"      ,   & ! hth
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         CALL histdef( nid_T, "so20chgt", "Depth of 20C isotherm"              , "m"      ,   & ! hd20
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         CALL histdef( nid_T, "so28chgt", "Depth of 28C isotherm"              , "m"      ,   & ! hd28
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         CALL histdef( nid_T, "sohtc300", "Heat content 300 m"                 , "W"      ,   & ! htc3
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
#endif

         IF( ln_cpl .AND. nn_ice == 2 ) THEN
         CALL histdef( nid_T,"soicetem" , "Ice Surface Temperature"            , "K"      ,   &  ! tn_ice
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         CALL histdef( nid_T,"soicealb" , "Ice Albedo"                         , "[0,1]"  ,   &  ! alb_ice
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         ENDIF

         CALL histend( nid_T, snc4chunks=snc4set )

         !                                                                                      !!! nid_U : 3D
         CALL histdef( nid_U, "vozocrtx", "Zonal Current"                      , "m/s"    ,   &  ! un
            &          jpi, jpj, nh_U, ipk, 1, ipk, nz_U, 32, clop, zsto, zout )
         IF( ln_traldf_gdia ) THEN
            CALL histdef( nid_U, "vozoeivu", "Zonal EIV Current"                  , "m/s"    ,   &  ! u_eiv
                 &          jpi, jpj, nh_U, ipk, 1, ipk, nz_U, 32, clop, zsto, zout )
         ELSE
#if defined key_diaeiv
         CALL histdef( nid_U, "vozoeivu", "Zonal EIV Current"                  , "m/s"    ,   &  ! u_eiv
            &          jpi, jpj, nh_U, ipk, 1, ipk, nz_U, 32, clop, zsto, zout )
#endif
         END IF
         !                                                                                      !!! nid_U : 2D
         CALL histdef( nid_U, "sozotaux", "Wind Stress along i-axis"           , "N/m2"   ,   &  ! utau
            &          jpi, jpj, nh_U, 1  , 1, 1  , - 99, 32, clop, zsto, zout )

         CALL histend( nid_U, snc4chunks=snc4set )

         !                                                                                      !!! nid_V : 3D
         CALL histdef( nid_V, "vomecrty", "Meridional Current"                 , "m/s"    ,   &  ! vn
            &          jpi, jpj, nh_V, ipk, 1, ipk, nz_V, 32, clop, zsto, zout )
         IF( ln_traldf_gdia ) THEN
            CALL histdef( nid_V, "vomeeivv", "Meridional EIV Current"             , "m/s"    ,   &  ! v_eiv
                 &          jpi, jpj, nh_V, ipk, 1, ipk, nz_V, 32, clop, zsto, zout )
         ELSE 
#if defined key_diaeiv
         CALL histdef( nid_V, "vomeeivv", "Meridional EIV Current"             , "m/s"    ,   &  ! v_eiv
            &          jpi, jpj, nh_V, ipk, 1, ipk, nz_V, 32, clop, zsto, zout )
#endif
         END IF
         !                                                                                      !!! nid_V : 2D
         CALL histdef( nid_V, "sometauy", "Wind Stress along j-axis"           , "N/m2"   ,   &  ! vtau
            &          jpi, jpj, nh_V, 1  , 1, 1  , - 99, 32, clop, zsto, zout )

         CALL histend( nid_V, snc4chunks=snc4set )

         !                                                                                      !!! nid_W : 3D
         CALL histdef( nid_W, "vovecrtz", "Vertical Velocity"                  , "m/s"    ,   &  ! wn
            &          jpi, jpj, nh_W, ipk, 1, ipk, nz_W, 32, clop, zsto, zout )
         IF( ln_traldf_gdia ) THEN
            CALL histdef( nid_W, "voveeivw", "Vertical EIV Velocity"              , "m/s"    ,   &  ! w_eiv
                 &          jpi, jpj, nh_W, ipk, 1, ipk, nz_W, 32, clop, zsto, zout )
         ELSE
#if defined key_diaeiv
         CALL histdef( nid_W, "voveeivw", "Vertical EIV Velocity"              , "m/s"    ,   &  ! w_eiv
            &          jpi, jpj, nh_W, ipk, 1, ipk, nz_W, 32, clop, zsto, zout )
#endif
         END IF
         CALL histdef( nid_W, "votkeavt", "Vertical Eddy Diffusivity"          , "m2/s"   ,   &  ! avt
            &          jpi, jpj, nh_W, ipk, 1, ipk, nz_W, 32, clop, zsto, zout )
         CALL histdef( nid_W, "votkeavm", "Vertical Eddy Viscosity"             , "m2/s"  ,   &  ! avmu
            &          jpi, jpj, nh_W, ipk, 1, ipk, nz_W, 32, clop, zsto, zout )

         IF( lk_zdfddm ) THEN
            CALL histdef( nid_W,"voddmavs","Salt Vertical Eddy Diffusivity"    , "m2/s"   ,   &  ! avs
               &          jpi, jpj, nh_W, ipk, 1, ipk, nz_W, 32, clop, zsto, zout )
         ENDIF
         !                                                                                      !!! nid_W : 2D
#if defined key_traldf_c2d
         CALL histdef( nid_W, "soleahtw", "lateral eddy diffusivity"           , "m2/s"   ,   &  ! ahtw
            &          jpi, jpj, nh_W, 1  , 1, 1  , - 99, 32, clop, zsto, zout )
# if defined key_traldf_eiv 
            CALL histdef( nid_W, "soleaeiw", "eddy induced vel. coeff. at w-point", "m2/s",   &  ! aeiw
               &       jpi, jpj, nh_W, 1  , 1, 1  , - 99, 32, clop, zsto, zout )
# endif
#endif

         CALL histend( nid_W, snc4chunks=snc4set )

         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'End of NetCDF Initialization'
         IF(ll_print) CALL FLUSH(numout )

      ENDIF

      ! 2. Start writing data
      ! ---------------------

      ! ndex(1) est utilise ssi l'avant dernier argument est different de 
      ! la taille du tableau en sortie. Dans ce cas , l'avant dernier argument
      ! donne le nombre d'elements, et ndex la liste des indices a sortir

      IF( lwp .AND. MOD( itmod, nwrite ) == 0 ) THEN 
         WRITE(numout,*) 'dia_wri : write model outputs in NetCDF files at ', kt, 'time-step'
         WRITE(numout,*) '~~~~~~ '
      ENDIF

      IF( lk_vvl ) THEN
         CALL histwrite( nid_T, "votemper", it, tsn(:,:,:,jp_tem) * fse3t_n(:,:,:) , ndim_T , ndex_T  )   ! heat content
         CALL histwrite( nid_T, "vosaline", it, tsn(:,:,:,jp_sal) * fse3t_n(:,:,:) , ndim_T , ndex_T  )   ! salt content
         CALL histwrite( nid_T, "sosstsst", it, tsn(:,:,1,jp_tem) * fse3t_n(:,:,1) , ndim_hT, ndex_hT )   ! sea surface heat content
         CALL histwrite( nid_T, "sosaline", it, tsn(:,:,1,jp_sal) * fse3t_n(:,:,1) , ndim_hT, ndex_hT )   ! sea surface salinity content
      ELSE
      CALL histwrite( nid_T, "votemper", it, tsn(:,:,:,jp_tem), ndim_T , ndex_T  )   ! temperature
      CALL histwrite( nid_T, "vosaline", it, tsn(:,:,:,jp_sal), ndim_T , ndex_T  )   ! salinity
      CALL histwrite( nid_T, "sosstsst", it, tsn(:,:,1,jp_tem), ndim_hT, ndex_hT )   ! sea surface temperature
      CALL histwrite( nid_T, "sosaline", it, tsn(:,:,1,jp_sal), ndim_hT, ndex_hT )   ! sea surface salinity
      ENDIF
      IF( lk_vvl ) THEN
         zw3d(:,:,:) = ( ( fse3t_n(:,:,:) - e3t_0(:,:,:) ) / e3t_0(:,:,:) * 100 * tmask(:,:,:) ) ** 2
         CALL histwrite( nid_T, "vovvle3t", it, fse3t_n (:,:,:) , ndim_T , ndex_T  )   ! level thickness
         CALL histwrite( nid_T, "vovvldep", it, fsdept_n(:,:,:) , ndim_T , ndex_T  )   ! t-point depth
         CALL histwrite( nid_T, "vovvldef", it, zw3d             , ndim_T , ndex_T  )   ! level thickness deformation
      ENDIF
      CALL histwrite( nid_T, "sossheig", it, sshn          , ndim_hT, ndex_hT )   ! sea surface height
      CALL histwrite( nid_T, "sowaflup", it, ( emp-rnf )   , ndim_hT, ndex_hT )   ! upward water flux
      CALL histwrite( nid_T, "sorunoff", it, rnf           , ndim_hT, ndex_hT )   ! river runoffs
      CALL histwrite( nid_T, "sosfldow", it, sfx           , ndim_hT, ndex_hT )   ! downward salt flux 
                                                                                  ! (includes virtual salt flux beneath ice 
                                                                                  ! in linear free surface case)
      IF( .NOT. lk_vvl ) THEN
      zw2d(:,:) = emp (:,:) * tsn(:,:,1,jp_tem)
      CALL histwrite( nid_T, "sosst_cd", it, zw2d, ndim_hT, ndex_hT )             ! c/d term on sst
      zw2d(:,:) = emp (:,:) * tsn(:,:,1,jp_sal)
      CALL histwrite( nid_T, "sosss_cd", it, zw2d, ndim_hT, ndex_hT )             ! c/d term on sss
      ENDIF
      CALL histwrite( nid_T, "sohefldo", it, qns + qsr     , ndim_hT, ndex_hT )   ! total heat flux
      CALL histwrite( nid_T, "soshfldo", it, qsr           , ndim_hT, ndex_hT )   ! solar heat flux
      CALL histwrite( nid_T, "somixhgt", it, hmld          , ndim_hT, ndex_hT )   ! turbocline depth
      CALL histwrite( nid_T, "somxl010", it, hmlp          , ndim_hT, ndex_hT )   ! mixed layer depth
      CALL histwrite( nid_T, "soicecov", it, fr_i          , ndim_hT, ndex_hT )   ! ice fraction   
      CALL histwrite( nid_T, "sowindsp", it, wndm          , ndim_hT, ndex_hT )   ! wind speed   
!
      IF( ln_icebergs ) THEN
         !
         CALL histwrite( nid_T, "calving"             , it, berg_grid%calving      , ndim_hT, ndex_hT )  
         CALL histwrite( nid_T, "calving_heat"        , it, berg_grid%calving_hflx , ndim_hT, ndex_hT )         
         CALL histwrite( nid_T, "berg_floating_melt"  , it, berg_grid%floating_melt, ndim_hT, ndex_hT )  
         !
         CALL histwrite( nid_T, "berg_stored_ice"     , it, berg_grid%stored_ice   , ndim_bT, ndex_bT )
         !
         IF( ln_bergdia ) THEN
            CALL histwrite( nid_T, "berg_melt"           , it, berg_melt        , ndim_hT, ndex_hT   )  
            CALL histwrite( nid_T, "berg_buoy_melt"      , it, buoy_melt        , ndim_hT, ndex_hT   )  
            CALL histwrite( nid_T, "berg_eros_melt"      , it, eros_melt        , ndim_hT, ndex_hT   )  
            CALL histwrite( nid_T, "berg_conv_melt"      , it, conv_melt        , ndim_hT, ndex_hT   )  
            CALL histwrite( nid_T, "berg_virtual_area"   , it, virtual_area     , ndim_hT, ndex_hT   )  
            CALL histwrite( nid_T, "bits_src"            , it, bits_src         , ndim_hT, ndex_hT   )  
            CALL histwrite( nid_T, "bits_melt"           , it, bits_melt        , ndim_hT, ndex_hT   )  
            CALL histwrite( nid_T, "bits_mass"           , it, bits_mass        , ndim_hT, ndex_hT   )  
            CALL histwrite( nid_T, "berg_mass"           , it, berg_mass        , ndim_hT, ndex_hT   )  
            !
            CALL histwrite( nid_T, "berg_real_calving"   , it, real_calving     , ndim_bT, ndex_bT   )
         ENDIF
      ENDIF

      IF( .NOT. ln_cpl ) THEN
      CALL histwrite( nid_T, "sohefldp", it, qrp           , ndim_hT, ndex_hT )   ! heat flux damping
      CALL histwrite( nid_T, "sowafldp", it, erp           , ndim_hT, ndex_hT )   ! freshwater flux damping
      IF( ln_ssr ) zw2d(:,:) = erp(:,:) * tsn(:,:,1,jp_sal) * tmask(:,:,1)
      CALL histwrite( nid_T, "sosafldp", it, zw2d          , ndim_hT, ndex_hT )   ! salt flux damping
      ENDIF
      IF( ln_cpl .AND. nn_ice <= 1 ) THEN
      CALL histwrite( nid_T, "sohefldp", it, qrp           , ndim_hT, ndex_hT )   ! heat flux damping
      CALL histwrite( nid_T, "sowafldp", it, erp           , ndim_hT, ndex_hT )   ! freshwater flux damping
      IF( ln_ssr ) zw2d(:,:) = erp(:,:) * tsn(:,:,1,jp_sal) * tmask(:,:,1)
      CALL histwrite( nid_T, "sosafldp", it, zw2d          , ndim_hT, ndex_hT )   ! salt flux damping
      ENDIF
!      zw2d(:,:) = FLOAT( nmln(:,:) ) * tmask(:,:,1)
!      CALL histwrite( nid_T, "sobowlin", it, zw2d          , ndim_hT, ndex_hT )   ! ???

#if defined key_diahth
      CALL histwrite( nid_T, "sothedep", it, hth           , ndim_hT, ndex_hT )   ! depth of the thermocline
      CALL histwrite( nid_T, "so20chgt", it, hd20          , ndim_hT, ndex_hT )   ! depth of the 20 isotherm
      CALL histwrite( nid_T, "so28chgt", it, hd28          , ndim_hT, ndex_hT )   ! depth of the 28 isotherm
      CALL histwrite( nid_T, "sohtc300", it, htc3          , ndim_hT, ndex_hT )   ! first 300m heaat content
#endif

      IF( ln_cpl .AND. nn_ice == 2 ) THEN
      CALL histwrite( nid_T, "soicetem", it, tn_ice(:,:,1) , ndim_hT, ndex_hT )   ! surf. ice temperature
      CALL histwrite( nid_T, "soicealb", it, alb_ice(:,:,1), ndim_hT, ndex_hT )   ! ice albedo
      ENDIF

      CALL histwrite( nid_U, "vozocrtx", it, un            , ndim_U , ndex_U )    ! i-current
      IF( ln_traldf_gdia ) THEN
         IF (.not. ALLOCATED(psix_eiv))THEN
            ALLOCATE( psix_eiv(jpi,jpj,jpk) , psiy_eiv(jpi,jpj,jpk) , STAT=ierr )
            IF( lk_mpp   )   CALL mpp_sum ( ierr )
            IF( ierr > 0 )   CALL ctl_stop('STOP', 'diawri: unable to allocate psi{x,y}_eiv')
            psix_eiv(:,:,:) = 0.0_wp
            psiy_eiv(:,:,:) = 0.0_wp
         ENDIF
         DO jk=1,jpkm1
            zw3d(:,:,jk) = (psix_eiv(:,:,jk+1) - psix_eiv(:,:,jk))/fse3u(:,:,jk)  ! u_eiv = -dpsix/dz
         END DO
         zw3d(:,:,jpk) = 0._wp
         CALL histwrite( nid_U, "vozoeivu", it, zw3d, ndim_U , ndex_U )           ! i-eiv current
      ELSE
#if defined key_diaeiv
      CALL histwrite( nid_U, "vozoeivu", it, u_eiv         , ndim_U , ndex_U )    ! i-eiv current
#endif
      ENDIF
      CALL histwrite( nid_U, "sozotaux", it, utau          , ndim_hU, ndex_hU )   ! i-wind stress

      CALL histwrite( nid_V, "vomecrty", it, vn            , ndim_V , ndex_V  )   ! j-current
      IF( ln_traldf_gdia ) THEN
         DO jk=1,jpk-1
            zw3d(:,:,jk) = (psiy_eiv(:,:,jk+1) - psiy_eiv(:,:,jk))/fse3v(:,:,jk)  ! v_eiv = -dpsiy/dz
         END DO
         zw3d(:,:,jpk) = 0._wp
         CALL histwrite( nid_V, "vomeeivv", it, zw3d, ndim_V , ndex_V )           ! j-eiv current
      ELSE
#if defined key_diaeiv
      CALL histwrite( nid_V, "vomeeivv", it, v_eiv         , ndim_V , ndex_V  )   ! j-eiv current
#endif
      ENDIF
      CALL histwrite( nid_V, "sometauy", it, vtau          , ndim_hV, ndex_hV )   ! j-wind stress

      CALL histwrite( nid_W, "vovecrtz", it, wn             , ndim_T, ndex_T )    ! vert. current
      IF( ln_traldf_gdia ) THEN
         DO jk=1,jpk-1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1  ! vector opt.
                  zw3d(ji,jj,jk) = (psiy_eiv(ji,jj,jk) - psiy_eiv(ji,jj-1,jk))/e2v(ji,jj) + &
                       &    (psix_eiv(ji,jj,jk) - psix_eiv(ji-1,jj,jk))/e1u(ji,jj) ! w_eiv = dpsiy/dy + dpsiy/dx
               END DO
            END DO
         END DO
         zw3d(:,:,jpk) = 0._wp
         CALL histwrite( nid_W, "voveeivw", it, zw3d          , ndim_T, ndex_T )    ! vert. eiv current
      ELSE
#   if defined key_diaeiv
      CALL histwrite( nid_W, "voveeivw", it, w_eiv          , ndim_T, ndex_T )    ! vert. eiv current
#   endif
      ENDIF
      CALL histwrite( nid_W, "votkeavt", it, avt            , ndim_T, ndex_T )    ! T vert. eddy diff. coef.
      CALL histwrite( nid_W, "votkeavm", it, avmu           , ndim_T, ndex_T )    ! T vert. eddy visc. coef.
      IF( lk_zdfddm ) THEN
         CALL histwrite( nid_W, "voddmavs", it, fsavs(:,:,:), ndim_T, ndex_T )    ! S vert. eddy diff. coef.
      ENDIF
#if defined key_traldf_c2d
      CALL histwrite( nid_W, "soleahtw", it, ahtw          , ndim_hT, ndex_hT )   ! lateral eddy diff. coef.
# if defined key_traldf_eiv
         CALL histwrite( nid_W, "soleaeiw", it, aeiw       , ndim_hT, ndex_hT )   ! EIV coefficient at w-point
# endif
#endif

      ! 3. Close all files
      ! ---------------------------------------
      IF( kt == nitend ) THEN
         CALL histclo( nid_T )
         CALL histclo( nid_U )
         CALL histclo( nid_V )
         CALL histclo( nid_W )
      ENDIF
      !
      CALL wrk_dealloc( jpi , jpj      , zw2d )
      IF ( ln_traldf_gdia .OR. lk_vvl )  call wrk_dealloc( jpi , jpj , jpk  , zw3d )
      !
      IF( nn_timing == 1 )   CALL timing_stop('dia_wri')
      !
   END SUBROUTINE dia_wri
# endif



   SUBROUTINE dia_wri_state( cdfile_name, kt )
      !!---------------------------------------------------------------------
      !!                 ***  ROUTINE dia_wri_state  ***
      !!        
      !! ** Purpose :   create a NetCDF file named cdfile_name which contains 
      !!      the instantaneous ocean state and forcing fields.
      !!        Used to find errors in the initial state or save the last
      !!      ocean state in case of abnormal end of a simulation
      !!
      !! ** Method  :   NetCDF files using ioipsl
      !!      File 'output.init.nc'  is created if ninist = 1 (namelist)
      !!      File 'output.abort.nc' is created in case of abnormal job end
      !!----------------------------------------------------------------------
      CHARACTER (len=* ), INTENT( in ) ::   cdfile_name      ! name of the file created
      INTEGER           , INTENT( in ) ::   kt               ! ocean time-step index
      !! 
      CHARACTER (len=32) :: clname
      CHARACTER (len=40) :: clop
      INTEGER  ::   id_i , nz_i, nh_i       
      INTEGER, DIMENSION(1) ::   idex             ! local workspace
      REAL(wp) ::   zsto, zout, zmax, zjulian, zdt
      !!----------------------------------------------------------------------
      ! 
!     IF( nn_timing == 1 )   CALL timing_start('dia_wri_state') ! not sure this works for routines not called in first timestep

      ! 0. Initialisation
      ! -----------------

      ! Define name, frequency of output and means
      clname = cdfile_name
      IF( .NOT. Agrif_Root() ) clname = TRIM(Agrif_CFixed())//'_'//TRIM(clname)
      zdt  = rdt
      zsto = rdt
      clop = "inst(x)"           ! no use of the mask value (require less cpu time)
      zout = rdt
      zmax = ( nitend - nit000 + 1 ) * zdt

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'dia_wri_state : single instantaneous ocean state'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~   and forcing fields file created '
      IF(lwp) WRITE(numout,*) '                and named :', clname, '.nc'


      ! 1. Define NETCDF files and fields at beginning of first time step
      ! -----------------------------------------------------------------

      ! Compute julian date from starting date of the run
      CALL ymds2ju( nyear, nmonth, nday, rdt, zjulian )         ! time axis 
      zjulian = zjulian - adatrj   !   set calendar origin to the beginning of the experiment
      CALL histbeg( clname, jpi, glamt, jpj, gphit,   &
          1, jpi, 1, jpj, nit000-1, zjulian, zdt, nh_i, id_i, domain_id=nidom, snc4chunks=snc4set ) ! Horizontal grid : glamt and gphit
      CALL histvert( id_i, "deptht", "Vertical T levels",   &    ! Vertical grid : gdept
          "m", jpk, gdept_1d, nz_i, "down")

      ! Declare all the output fields as NetCDF variables

      CALL histdef( id_i, "vosaline", "Salinity"              , "PSU"    ,   &   ! salinity
         &          jpi, jpj, nh_i, jpk, 1, jpk, nz_i, 32, clop, zsto, zout )
      CALL histdef( id_i, "votemper", "Temperature"           , "C"      ,   &   ! temperature
         &          jpi, jpj, nh_i, jpk, 1, jpk, nz_i, 32, clop, zsto, zout )
      CALL histdef( id_i, "sossheig", "Sea Surface Height"    , "m"      ,   &  ! ssh
         &          jpi, jpj, nh_i, 1  , 1, 1  , nz_i, 32, clop, zsto, zout )
      CALL histdef( id_i, "vozocrtx", "Zonal Current"         , "m/s"    ,   &   ! zonal current
         &          jpi, jpj, nh_i, jpk, 1, jpk, nz_i, 32, clop, zsto, zout )
      CALL histdef( id_i, "vomecrty", "Meridional Current"    , "m/s"    ,   &   ! meridonal current
         &          jpi, jpj, nh_i, jpk, 1, jpk, nz_i, 32, clop, zsto, zout ) 
      CALL histdef( id_i, "vovecrtz", "Vertical Velocity"     , "m/s"    ,   &   ! vertical current
         &          jpi, jpj, nh_i, jpk, 1, jpk, nz_i, 32, clop, zsto, zout ) 
      CALL histdef( id_i, "sowaflup", "Net Upward Water Flux" , "Kg/m2/S",   &   ! net freshwater 
         &          jpi, jpj, nh_i, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
      CALL histdef( id_i, "sohefldo", "Net Downward Heat Flux", "W/m2"   ,   &   ! net heat flux
         &          jpi, jpj, nh_i, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
      CALL histdef( id_i, "soshfldo", "Shortwave Radiation"   , "W/m2"   ,   &   ! solar flux
         &          jpi, jpj, nh_i, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
      CALL histdef( id_i, "soicecov", "Ice fraction"          , "[0,1]"  ,   &   ! fr_i
         &          jpi, jpj, nh_i, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
      CALL histdef( id_i, "sozotaux", "Zonal Wind Stress"     , "N/m2"   ,   &   ! i-wind stress
         &          jpi, jpj, nh_i, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
      CALL histdef( id_i, "sometauy", "Meridional Wind Stress", "N/m2"   ,   &   ! j-wind stress
         &          jpi, jpj, nh_i, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
      IF( lk_vvl ) THEN
         CALL histdef( id_i, "vovvldep", "T point depth"         , "m"      ,   &   ! t-point depth
            &          jpi, jpj, nh_i, jpk, 1, jpk, nz_i, 32, clop, zsto, zout )
         CALL histdef( id_i, "vovvle3t", "T point thickness"         , "m"      ,   &   ! t-point depth
            &          jpi, jpj, nh_i, jpk, 1, jpk, nz_i, 32, clop, zsto, zout )
      END IF

#if defined key_lim2
      CALL lim_wri_state_2( kt, id_i, nh_i )
#elif defined key_lim3
      CALL lim_wri_state( kt, id_i, nh_i )
#else
      CALL histend( id_i, snc4chunks=snc4set )
#endif

      ! 2. Start writing data
      ! ---------------------
      ! idex(1) est utilise ssi l'avant dernier argument est diffferent de 
      ! la taille du tableau en sortie. Dans ce cas , l'avant dernier argument
      ! donne le nombre d'elements, et idex la liste des indices a sortir
      idex(1) = 1   ! init to avoid compil warning

      ! Write all fields on T grid
      CALL histwrite( id_i, "votemper", kt, tsn(:,:,:,jp_tem), jpi*jpj*jpk, idex )    ! now temperature
      CALL histwrite( id_i, "vosaline", kt, tsn(:,:,:,jp_sal), jpi*jpj*jpk, idex )    ! now salinity
      CALL histwrite( id_i, "sossheig", kt, sshn             , jpi*jpj    , idex )    ! sea surface height
      CALL histwrite( id_i, "vozocrtx", kt, un               , jpi*jpj*jpk, idex )    ! now i-velocity
      CALL histwrite( id_i, "vomecrty", kt, vn               , jpi*jpj*jpk, idex )    ! now j-velocity
      CALL histwrite( id_i, "vovecrtz", kt, wn               , jpi*jpj*jpk, idex )    ! now k-velocity
      CALL histwrite( id_i, "sowaflup", kt, (emp-rnf )       , jpi*jpj    , idex )    ! freshwater budget
      CALL histwrite( id_i, "sohefldo", kt, qsr + qns        , jpi*jpj    , idex )    ! total heat flux
      CALL histwrite( id_i, "soshfldo", kt, qsr              , jpi*jpj    , idex )    ! solar heat flux
      CALL histwrite( id_i, "soicecov", kt, fr_i             , jpi*jpj    , idex )    ! ice fraction
      CALL histwrite( id_i, "sozotaux", kt, utau             , jpi*jpj    , idex )    ! i-wind stress
      CALL histwrite( id_i, "sometauy", kt, vtau             , jpi*jpj    , idex )    ! j-wind stress
      IF( lk_vvl ) THEN
         CALL histwrite( id_i, "vovvldep", kt, fsdept_n(:,:,:), jpi*jpj*jpk, idex )!  T-cell depth       
         CALL histwrite( id_i, "vovvle3t", kt, fse3t_n (:,:,:), jpi*jpj*jpk, idex )!  T-cell thickness  
      END IF

      ! 3. Close the file
      ! -----------------
      CALL histclo( id_i )
#if ! defined key_iomput && ! defined key_dimgout
      IF( ninist /= 1  ) THEN
         CALL histclo( nid_T )
         CALL histclo( nid_U )
         CALL histclo( nid_V )
         CALL histclo( nid_W )
      ENDIF
#endif
       
!     IF( nn_timing == 1 )   CALL timing_stop('dia_wri_state') ! not sure this works for routines not called in first timestep
      ! 

   END SUBROUTINE dia_wri_state
#endif
   !!======================================================================
END MODULE diawri
