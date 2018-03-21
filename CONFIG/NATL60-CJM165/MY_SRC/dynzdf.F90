MODULE dynzdf
   !!==============================================================================
   !!                 ***  MODULE  dynzdf  ***
   !! Ocean dynamics :  vertical component of the momentum mixing trend
   !!==============================================================================
   !! History :  1.0  !  2005-11  (G. Madec)  Original code
   !!            3.3  !  2010-10  (C. Ethe, G. Madec) reorganisation of initialisation phase
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dyn_zdf      : Update the momentum trend with the vertical diffusion
   !!   dyn_zdf_init : initializations of the vertical diffusion scheme
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain variables 
   USE zdf_oce         ! ocean vertical physics variables
   USE iom            ! I/O manager library

   USE dynzdf_exp      ! vertical diffusion: explicit (dyn_zdf_exp     routine)
   USE dynzdf_imp      ! vertical diffusion: implicit (dyn_zdf_imp     routine)

   USE ldfdyn_oce      ! ocean dynamics: lateral physics
   USE trd_oce         ! trends: ocean variables
   USE trddyn          ! trend manager: dynamics
   USE diavecq         ! diags for Q vector
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! MPP library
   USE prtctl          ! Print control
   USE wrk_nemo        ! Memory Allocation
   USE timing          ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dyn_zdf       !  routine called by step.F90
   PUBLIC   dyn_zdf_init  !  routine called by opa.F90

   INTEGER  ::   nzdf = 0   ! type vertical diffusion algorithm used, defined from ln_zdf... namlist logicals
   REAL(wp) ::   r2dt       ! time-step, = 2 rdttra except at nit000 (=rdttra) if neuler=0

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "zdfddm_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: dynzdf.F90 4990 2014-12-15 16:42:49Z timgraham $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS
   
   SUBROUTINE dyn_zdf( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_zdf  ***
      !!
      !! ** Purpose :   compute the vertical ocean dynamics physics.
      !!---------------------------------------------------------------------
      !!
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
      !
      INTEGER  ::   ji, jj, jk                      ! dummy loop arguments
      REAL(wp), POINTER, DIMENSION(:,:,:) ::  ztrdu, ztrdv
      REAL(wp), POINTER, DIMENSION(:,:,:) ::  zdzuzdf, zdzvzdf
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('dyn_zdf')
      !
      !                                          ! set time step
      IF( neuler == 0 .AND. kt == nit000     ) THEN   ;   r2dt =      rdt   ! = rdtra (restart with Euler time stepping)
      ELSEIF(               kt <= nit000 + 1 ) THEN   ;   r2dt = 2. * rdt   ! = 2 rdttra (leapfrog)
      ENDIF

      IF( l_trddyn )   THEN                      ! temporary save of ta and sa trends
         CALL wrk_alloc( jpi, jpj, jpk, ztrdu, ztrdv ) 
         ztrdu(:,:,:) = ua(:,:,:)
         ztrdv(:,:,:) = va(:,:,:)
      ENDIF
      IF( lk_diavecq )   THEN                      ! temporary save of ta and sa trends
         CALL wrk_alloc( jpi, jpj, jpk, ztrdu, ztrdv )
         CALL wrk_alloc( jpi, jpj, jpk, zdzuzdf, zdzvzdf )
         ztrdu(:,:,:) = ua(:,:,:)
         ztrdv(:,:,:) = va(:,:,:)
         zdzuzdf(:,:,:) = 0.
         zdzvzdf(:,:,:) = 0.
      ENDIF

      SELECT CASE ( nzdf )                       ! compute lateral mixing trend and add it to the general trend
      !
      CASE ( 0 )   ;   CALL dyn_zdf_exp( kt, r2dt )      ! explicit scheme
      CASE ( 1 )   ;   CALL dyn_zdf_imp( kt, r2dt )      ! implicit scheme
      !
      CASE ( -1 )                                        ! esopa: test all possibility with control print
                       CALL dyn_zdf_exp( kt, r2dt )
                       CALL prt_ctl( tab3d_1=ua, clinfo1=' zdf0 - Ua: ', mask1=umask,               &
                          &          tab3d_2=va, clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
                       CALL dyn_zdf_imp( kt, r2dt )
                       CALL prt_ctl( tab3d_1=ua, clinfo1=' zdf1 - Ua: ', mask1=umask,               &
                          &          tab3d_2=va, clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
      END SELECT

      IF( l_trddyn )   THEN                      ! save the vertical diffusive trends for further diagnostics
         ztrdu(:,:,:) = ua(:,:,:) - ztrdu(:,:,:)
         ztrdv(:,:,:) = va(:,:,:) - ztrdv(:,:,:)
         CALL trd_dyn( ztrdu, ztrdv, jpdyn_zdf, kt )
         CALL wrk_dealloc( jpi, jpj, jpk, ztrdu, ztrdv ) 
      ENDIF
      !
      IF( lk_diavecq )   THEN                     ! save the vertical flux of momentum due to viscosity
         ztrdu(:,:,:) = ua(:,:,:) - ztrdu(:,:,:)
         ztrdv(:,:,:) = va(:,:,:) - ztrdv(:,:,:)

         DO jk = 1, jpkm1
           DO jj = 1, jpjm1
             DO ji = 1, jpim1 

               zdzuzdf(ji,jj,jk) = ( ztrdu(ji,jj,jk+1) - ztrdu(ji,jj,jk) ) * 2. / ( e3w_0(ji+1,jj,jk) + e3w_0(ji,jj,jk) )
               zdzvzdf(ji,jj,jk) = ( ztrdv(ji,jj,jk+1) - ztrdv(ji,jj,jk) ) * 2. / ( e3w_0(ji,jj+1,jk) + e3w_0(ji,jj,jk) )

             END DO
           END DO
         END DO

         CALL iom_put( 'dzuzdf', zdzuzdf  ) 
         CALL iom_put( 'dzvzdf', zdzvzdf  ) 
         CALL wrk_dealloc( jpi, jpj, jpk, ztrdu, ztrdv )
         CALL wrk_dealloc( jpi, jpj, jpk, zdzuzdf, zdzvzdf )
      ENDIF
      !                                          ! print mean trends (used for debugging)
      IF(ln_ctl)   CALL prt_ctl( tab3d_1=ua, clinfo1=' zdf  - Ua: ', mask1=umask,               &
            &                    tab3d_2=va, clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
      !
      IF( nn_timing == 1 )  CALL timing_stop('dyn_zdf')
      !
   END SUBROUTINE dyn_zdf


   SUBROUTINE dyn_zdf_init
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE dyn_zdf_init  ***
      !!
      !! ** Purpose :   initializations of the vertical diffusion scheme
      !!
      !! ** Method  :   implicit (euler backward) scheme (default)
      !!                explicit (time-splitting) scheme if ln_zdfexp=T
      !!----------------------------------------------------------------------
      USE zdftke
      USE zdfgls
      USE zdfkpp
      !!----------------------------------------------------------------------
      !
      ! Choice from ln_zdfexp read in namelist in zdfini
      IF( ln_zdfexp ) THEN   ;   nzdf = 0           ! use explicit scheme
      ELSE                   ;   nzdf = 1           ! use implicit scheme
      ENDIF
      !
      ! Force implicit schemes
      IF( lk_zdftke .OR. lk_zdfgls .OR. lk_zdfkpp )   nzdf = 1   ! TKE, GLS or KPP physics
      IF( ln_dynldf_iso                           )   nzdf = 1   ! iso-neutral lateral physics
      IF( ln_dynldf_hor .AND. ln_sco              )   nzdf = 1   ! horizontal lateral physics in s-coordinate
      !
      IF( lk_esopa )    nzdf = -1                   ! Esopa key: All schemes used
      !
      IF(lwp) THEN                                  ! Print the choice
         WRITE(numout,*)
         WRITE(numout,*) 'dyn_zdf_init : vertical dynamics physics scheme'
         WRITE(numout,*) '~~~~~~~~~~~'
         IF( nzdf == -1 )   WRITE(numout,*) '              ESOPA test All scheme used'
         IF( nzdf ==  0 )   WRITE(numout,*) '              Explicit time-splitting scheme'
         IF( nzdf ==  1 )   WRITE(numout,*) '              Implicit (euler backward) scheme'
      ENDIF
      !
   END SUBROUTINE dyn_zdf_init

   !!==============================================================================
END MODULE dynzdf
