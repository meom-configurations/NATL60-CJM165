MODULE trazdf
   !!==============================================================================
   !!                 ***  MODULE  trazdf  ***
   !! Ocean active tracers:  vertical component of the tracer mixing trend
   !!==============================================================================
   !! History :  1.0  ! 2005-11  (G. Madec)  Original code
   !!            3.0  ! 2008-01  (C. Ethe, G. Madec)  merge TRC-TRA
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   tra_zdf      : Update the tracer trend with the vertical diffusion
   !!   tra_zdf_init : initialisation of the computation
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain variables 
   USE domvvl          ! variable volume
   USE phycst          ! physical constant
   USE zdf_oce         ! ocean vertical physics variables
   USE sbc_oce         ! surface boundary condition: ocean
   USE dynspg_oce
   USE trazdf_exp      ! vertical diffusion: explicit (tra_zdf_exp     routine)
   USE trazdf_imp      ! vertical diffusion: implicit (tra_zdf_imp     routine)
   USE ldftra_oce      ! ocean active tracers: lateral physics
   USE trd_oce         ! trends: ocean variables
   USE trdtra          ! trends manager: tracers 
   USE diavecq         ! diags for Q vector
   USE eosbn2         ! equation of state                (eos_bn2 routine)
   USE iom            ! I/O manager library
   !
   USE in_out_manager  ! I/O manager
   USE prtctl          ! Print control
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp         ! MPP library
   USE wrk_nemo        ! Memory allocation
   USE timing          ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tra_zdf        ! routine called by step.F90
   PUBLIC   tra_zdf_init   ! routine called by nemogcm.F90

   INTEGER ::   nzdf = 0   ! type vertical diffusion algorithm used (defined from ln_zdf...  namlist logicals)

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "zdfddm_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.7 , NEMO Consortium (2014)
   !! $Id: trazdf.F90 5385 2015-06-09 13:50:42Z cetlod $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE tra_zdf( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_zdf  ***
      !!
      !! ** Purpose :   compute the vertical ocean tracer physics.
      !!---------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
      !!
      INTEGER  ::   ji,jj,jk                   ! Dummy loop indices
      REAL(wp) :: zstp
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   ztrdt, ztrds   ! 3D workspace
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   zdxtzdf, zdxszdf   ! 3D workspace
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   zdytzdf, zdyszdf   ! 3D workspace
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   zdxbzdf, zdybzdf   ! 3D workspace
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   zalpha, zbeta   ! 3D workspace
      REAL(wp), POINTER, DIMENSION(:,:,:,:) ::   zpab   ! 4D workspace
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   zrhd   ! 4D workspace
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('tra_zdf')
      zstp=kt-nit000+1
      !
      IF( neuler == 0 .AND. kt == nit000 ) THEN     ! at nit000
         r2dtra(:) =  rdttra(:)                          ! = rdtra (restarting with Euler time stepping)
      ELSEIF( kt <= nit000 + 1) THEN                ! at nit000 or nit000+1
         r2dtra(:) = 2. * rdttra(:)                      ! = 2 rdttra (leapfrog)
      ENDIF

      IF( l_trdtra )   THEN                    !* Save ta and sa trends
         CALL wrk_alloc( jpi, jpj, jpk, ztrdt, ztrds )
         ztrdt(:,:,:) = tsa(:,:,:,jp_tem)
         ztrds(:,:,:) = tsa(:,:,:,jp_sal)
      ENDIF

      IF( lk_diavecq )   THEN                    !* Save ta and sa trends
      IF ( l_diavecq_out)  THEN
         CALL wrk_alloc( jpi, jpj, jpk, ztrdt, ztrds )
         CALL wrk_alloc( jpi, jpj, jpk, zdxtzdf, zdxszdf )
         CALL wrk_alloc( jpi, jpj, jpk, zdytzdf, zdyszdf )
         CALL wrk_alloc( jpi, jpj, jpk, zdxbzdf, zdybzdf )
         CALL wrk_alloc( jpi, jpj, jpk,jpts, zpab )
         CALL wrk_alloc( jpi, jpj, jpk, zalpha, zbeta )
         CALL wrk_alloc( jpi, jpj, jpk, zrhd )

         ztrdt(:,:,:) = tsa(:,:,:,jp_tem)
         ztrds(:,:,:) = tsa(:,:,:,jp_sal)
         zdxtzdf(:,:,:) = 0.
         zdxszdf(:,:,:) = 0.
         zdytzdf(:,:,:) = 0.
         zdyszdf(:,:,:) = 0.
         zpab(:,:,:,:) = 0.
         zalpha(:,:,:) = 0.
         zbeta(:,:,:) = 0.
         zdxbzdf(:,:,:) = 0.
         zdybzdf(:,:,:) = 0.
         zrhd(:,:,:) = 0.
      ENDIF
      ENDIF

      SELECT CASE ( nzdf )                       ! compute lateral mixing trend and add it to the general trend
      CASE ( 0 )    ;    CALL tra_zdf_exp( kt, nit000, 'TRA', r2dtra, nn_zdfexp, tsb, tsa, jpts )  !   explicit scheme 
      CASE ( 1 )    ;    CALL tra_zdf_imp( kt, nit000, 'TRA', r2dtra,            tsb, tsa, jpts )  !   implicit scheme 
      CASE ( -1 )                                       ! esopa: test all possibility with control print
         CALL tra_zdf_exp( kt, nit000, 'TRA', r2dtra, nn_zdfexp, tsb, tsa, jpts )
         CALL prt_ctl( tab3d_1=tsa(:,:,:,jp_tem), clinfo1=' zdf0 - Ta: ', mask1=tmask,               &
         &             tab3d_2=tsa(:,:,:,jp_sal), clinfo2=       ' Sa: ', mask2=tmask, clinfo3='tra' )
         CALL tra_zdf_imp( kt, nit000, 'TRA', r2dtra,            tsb, tsa, jpts ) 
         CALL prt_ctl( tab3d_1=tsa(:,:,:,jp_tem), clinfo1=' zdf1 - Ta: ', mask1=tmask,               &
         &             tab3d_2=tsa(:,:,:,jp_sal), clinfo2=       ' Sa: ', mask2=tmask, clinfo3='tra' )
      END SELECT
      ! DRAKKAR SSS control {
      ! JMM avoid negative salinities near river outlet ! Ugly fix
      ! JMM : restore negative salinities to small salinities:
      WHERE ( tsa(:,:,:,jp_sal) < 0._wp )   tsa(:,:,:,jp_sal) = 0.1_wp

      IF( l_trdtra )   THEN                      ! save the vertical diffusive trends for further diagnostics
         DO jk = 1, jpkm1
            ztrdt(:,:,jk) = ( ( tsa(:,:,jk,jp_tem) - tsb(:,:,jk,jp_tem) ) / r2dtra(jk) ) - ztrdt(:,:,jk)
            ztrds(:,:,jk) = ( ( tsa(:,:,jk,jp_sal) - tsb(:,:,jk,jp_sal) ) / r2dtra(jk) ) - ztrds(:,:,jk)
         END DO
         CALL lbc_lnk( ztrdt, 'T', 1. )
         CALL lbc_lnk( ztrds, 'T', 1. )
         CALL trd_tra( kt, 'TRA', jp_tem, jptra_zdf, ztrdt )
         CALL trd_tra( kt, 'TRA', jp_sal, jptra_zdf, ztrds )
         CALL wrk_dealloc( jpi, jpj, jpk, ztrdt, ztrds )
      ENDIF
      IF( lk_diavecq )   THEN                      ! save the vertical diffusive trends for further diagnostics
      IF ( l_diavecq_out)  THEN
         DO jk = 1, jpkm1
            ztrdt(:,:,jk) = ( ( tsa(:,:,jk,jp_tem) - tsb(:,:,jk,jp_tem) ) / r2dtra(jk) ) - ztrdt(:,:,jk)
            ztrds(:,:,jk) = ( ( tsa(:,:,jk,jp_sal) - tsb(:,:,jk,jp_sal) ) / r2dtra(jk) ) - ztrds(:,:,jk)
         END DO
         CALL lbc_lnk( ztrdt, 'T', 1. )
         CALL lbc_lnk( ztrds, 'T', 1. )

         CALL eos_rab(tsa,zpab)  ! ideally, we could have used an average now/after 
         CALL eos( tsa, zrhd, fsdept_n(:,:,:) )

         zalpha(:,:,:) = zpab(:,:,:,jp_tem)
         zbeta (:,:,:) = zpab(:,:,:,jp_sal)

         DO jk = 1, jpk
           DO jj = 1, jpj
             DO ji = 2, jpim1
               ! TODO interpoler alpha / beta at u,v points
               ! TODO : ne pas mettre alpha et beta dans la derivee.
               zdxtzdf(ji,jj,jk) = ( zalpha(ji+1,jj,jk) + zalpha(ji,jj,jk) ) * &
                   &               (  ztrdt(ji+1,jj,jk) -  ztrdt(ji-1,jj,jk) ) / e1u(ji,jj) 
               zdxszdf(ji,jj,jk) = ( zbeta (ji+1,jj,jk) + zbeta (ji,jj,jk) ) * &
                   &               (  ztrds(ji+1,jj,jk) -  ztrds(ji-1,jj,jk) ) / e1u(ji,jj) 
               zdxbzdf(ji,jj,jk) = -1 * grav * zrhd(ji,jj,jk) * ( zdxszdf(ji,jj,jk) - zdxtzdf(ji,jj,jk) ) / rau0
             END DO
           END DO

           DO jj = 2, jpjm1
             DO ji = 1, jpi
               zdytzdf(ji,jj,jk)= ( zalpha(ji,jj+1,jk) + zalpha(ji,jj,jk) ) * &
                   &              (  ztrdt(ji,jj-1,jk) -  ztrdt(ji,jj-1,jk) ) / e2v(ji,jj) 
               zdyszdf(ji,jj,jk)= ( zbeta (ji,jj+1,jk) + zbeta (ji,jj,jk) ) * &
                   &              (  ztrds(ji,jj-1,jk) -  ztrds(ji,jj-1,jk) ) / e2v(ji,jj) 
               zdybzdf(ji,jj,jk) = -1 * grav * zrhd(ji,jj,jk) * ( zdyszdf(ji,jj,jk) - zdytzdf(ji,jj,jk) ) / rau0
             END DO
           END DO
         END DO

         CALL lbc_lnk( zdxbzdf, 'U', 1. )
         CALL lbc_lnk( zdybzdf, 'V', 1. )

         CALL iom_put( 'dxbzdf', zdxbzdf ) 
         CALL iom_put( 'dybzdf', zdybzdf ) 

         CALL wrk_dealloc( jpi, jpj, jpk, ztrdt, ztrds )
         CALL wrk_dealloc( jpi, jpj, jpk, zdxtzdf, zdxszdf )
         CALL wrk_dealloc( jpi, jpj, jpk, zdytzdf, zdyszdf )
         CALL wrk_dealloc( jpi, jpj, jpk, zdxbzdf, zdybzdf )
         CALL wrk_dealloc( jpi, jpj, jpk,jpts, zpab )
         CALL wrk_dealloc( jpi, jpj, jpk, zrhd )
         CALL wrk_dealloc( jpi, jpj, jpk, zalpha, zbeta )
      ENDIF
      ENDIF

      !                                          ! print mean trends (used for debugging)
      IF(ln_ctl)   CALL prt_ctl( tab3d_1=tsa(:,:,:,jp_tem), clinfo1=' zdf  - Ta: ', mask1=tmask,               &
         &                       tab3d_2=tsa(:,:,:,jp_sal), clinfo2=       ' Sa: ', mask2=tmask, clinfo3='tra' )
      !
      IF( nn_timing == 1 )  CALL timing_stop('tra_zdf')
      !
   END SUBROUTINE tra_zdf


   SUBROUTINE tra_zdf_init
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE tra_zdf_init  ***
      !!
      !! ** Purpose :   Choose the vertical mixing scheme
      !!
      !! ** Method  :   Set nzdf from ln_zdfexp
      !!      nzdf = 0   explicit (time-splitting) scheme (ln_zdfexp=T)
      !!           = 1   implicit (euler backward) scheme (ln_zdfexp=F)
      !!      NB: rotation of lateral mixing operator or TKE or KPP scheme,
      !!      the implicit scheme is required.
      !!----------------------------------------------------------------------
      USE zdftke
      USE zdfgls
      USE zdfkpp
      !!----------------------------------------------------------------------

      ! Choice from ln_zdfexp already read in namelist in zdfini module
      IF( ln_zdfexp ) THEN   ;   nzdf = 0           ! use explicit scheme
      ELSE                   ;   nzdf = 1           ! use implicit scheme
      ENDIF

      ! Force implicit schemes
      IF( lk_zdftke .OR. lk_zdfgls .OR. lk_zdfkpp )   nzdf = 1      ! TKE, GLS or KPP physics
      IF( ln_traldf_iso                           )   nzdf = 1      ! iso-neutral lateral physics
      IF( ln_traldf_hor .AND. ln_sco              )   nzdf = 1      ! horizontal lateral physics in s-coordinate
      IF( ln_zdfexp .AND. nzdf == 1 )   CALL ctl_stop( 'tra_zdf : If using the rotation of lateral mixing operator',   &
            &                         ' TKE or KPP scheme, the implicit scheme is required, set ln_zdfexp = .false.' )

      ! Test: esopa
      IF( lk_esopa )    nzdf = -1                      ! All schemes used

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'tra_zdf_init : vertical tracer physics scheme'
         WRITE(numout,*) '~~~~~~~~~~~'
         IF( nzdf == -1 )   WRITE(numout,*) '              ESOPA test All scheme used'
         IF( nzdf ==  0 )   WRITE(numout,*) '              Explicit time-splitting scheme'
         IF( nzdf ==  1 )   WRITE(numout,*) '              Implicit (euler backward) scheme'
      ENDIF
      !
   END SUBROUTINE tra_zdf_init

   !!==============================================================================
END MODULE trazdf
