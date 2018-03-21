MODULE dynnxt
   !!=========================================================================
   !!                       ***  MODULE  dynnxt  ***
   !! Ocean dynamics: time stepping
   !!=========================================================================
   !! History :  OPA  !  1987-02  (P. Andrich, D. L Hostis)  Original code
   !!                 !  1990-10  (C. Levy, G. Madec)
   !!            7.0  !  1993-03  (M. Guyon)  symetrical conditions
   !!            8.0  !  1997-02  (G. Madec & M. Imbard)  opa, release 8.0
   !!            8.2  !  1997-04  (A. Weaver)  Euler forward step
   !!             -   !  1997-06  (G. Madec)  lateral boudary cond., lbc routine
   !!    NEMO    1.0  !  2002-08  (G. Madec)  F90: Free form and module
   !!             -   !  2002-10  (C. Talandier, A-M. Treguier) Open boundary cond.
   !!            2.0  !  2005-11  (V. Garnier) Surface pressure gradient organization
   !!            2.3  !  2007-07  (D. Storkey) Calls to BDY routines. 
   !!            3.2  !  2009-06  (G. Madec, R.Benshila)  re-introduce the vvl option
   !!            3.3  !  2010-09  (D. Storkey, E.O'Dea) Bug fix for BDY module
   !!            3.3  !  2011-03  (P. Oddo) Bug fix for time-splitting+(BDY-OBC) and not VVL
   !!            3.5  !  2013-07  (J. Chanut) Compliant with time splitting changes
   !!            3.7  !  2014-04  (G. Madec) add the diagnostic of the time filter trends
   !!-------------------------------------------------------------------------
  
   !!-------------------------------------------------------------------------
   !!   dyn_nxt      : obtain the next (after) horizontal velocity
   !!-------------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE sbc_oce         ! Surface boundary condition: ocean fields
   USE phycst          ! physical constants
   USE dynspg_oce      ! type of surface pressure gradient
   USE dynadv          ! dynamics: vector invariant versus flux form
   USE domvvl          ! variable volume
   USE bdy_oce         ! ocean open boundary conditions
   USE bdydta          ! ocean open boundary conditions
   USE bdydyn          ! ocean open boundary conditions
   USE bdyvol          ! ocean open boundary condition (bdy_vol routines)
   USE trd_oce         ! trends: ocean variables
   USE diavecq         ! diags for Q vector
   USE trddyn          ! trend manager: dynamics
   USE trdken          ! trend manager: kinetic energy
   !
   USE in_out_manager  ! I/O manager
   USE iom             ! I/O manager library
   USE lbclnk          ! lateral boundary condition (or mpp link)
   USE lib_mpp         ! MPP library
   USE wrk_nemo        ! Memory Allocation
   USE prtctl          ! Print control
   USE timing          ! Timing
#if defined key_agrif
   USE agrif_opa_interp
#endif

   IMPLICIT NONE
   PRIVATE

   PUBLIC    dyn_nxt   ! routine called by step.F90

   !! * Substitutions
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: dynnxt.F90 5628 2015-07-22 20:26:35Z mathiot $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dyn_nxt ( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_nxt  ***
      !!                   
      !! ** Purpose :   Compute the after horizontal velocity. Apply the boundary 
      !!             condition on the after velocity, achieved the time stepping 
      !!             by applying the Asselin filter on now fields and swapping 
      !!             the fields.
      !!
      !! ** Method  : * After velocity is compute using a leap-frog scheme:
      !!                       (ua,va) = (ub,vb) + 2 rdt (ua,va)
      !!             Note that with flux form advection and variable volume layer
      !!             (lk_vvl=T), the leap-frog is applied on thickness weighted
      !!             velocity.
      !!             Note also that in filtered free surface (lk_dynspg_flt=T),
      !!             the time stepping has already been done in dynspg module
      !!
      !!              * Apply lateral boundary conditions on after velocity 
      !!             at the local domain boundaries through lbc_lnk call,
      !!             at the one-way open boundaries (lk_bdy=T),
      !!             at the AGRIF zoom   boundaries (lk_agrif=T)
      !!
      !!              * Apply the time filter applied and swap of the dynamics
      !!             arrays to start the next time step:
      !!                (ub,vb) = (un,vn) + atfp [ (ub,vb) + (ua,va) - 2 (un,vn) ]
      !!                (un,vn) = (ua,va).
      !!             Note that with flux form advection and variable volume layer
      !!             (lk_vvl=T), the time filter is applied on thickness weighted
      !!             velocity.
      !!
      !! ** Action :   ub,vb   filtered before horizontal velocity of next time-step
      !!               un,vn   now horizontal velocity of next time-step
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      INTEGER  ::   iku, ikv     ! local integers
#if ! defined key_dynspg_flt
      REAL(wp) ::   z2dt         ! temporary scalar
#endif
      REAL(wp) ::   zue3a, zue3n, zue3b, zuf, zec      ! local scalars
      REAL(wp) ::   zve3a, zve3n, zve3b, zvf, z1_2dt   !   -      -
      REAL(wp), POINTER, DIMENSION(:,:)   ::  zue, zve
      REAL(wp), POINTER, DIMENSION(:,:,:) ::  ze3u_f, ze3v_f, zua, zva 
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )   CALL timing_start('dyn_nxt')
      !
      CALL wrk_alloc( jpi,jpj,jpk,  ze3u_f, ze3v_f, zua, zva )
      IF( lk_dynspg_ts )   CALL wrk_alloc( jpi,jpj, zue, zve )
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_nxt : time stepping'
         IF(lwp) WRITE(numout,*) '~~~~~~~'
      ENDIF

#if defined key_dynspg_flt
      !
      ! Next velocity :   Leap-frog time stepping already done in dynspg_flt.F routine
      ! -------------

      ! Update after velocity on domain lateral boundaries      (only local domain required)
      ! --------------------------------------------------
      CALL lbc_lnk( ua, 'U', -1. )         ! local domain boundaries
      CALL lbc_lnk( va, 'V', -1. ) 
      !
#else

# if defined key_dynspg_exp
      ! Next velocity :   Leap-frog time stepping
      ! -------------
      z2dt = 2. * rdt                                 ! Euler or leap-frog time step 
      IF( neuler == 0 .AND. kt == nit000 )  z2dt = rdt
      !
      IF( ln_dynadv_vec .OR. .NOT. lk_vvl ) THEN      ! applied on velocity
         DO jk = 1, jpkm1
            ua(:,:,jk) = ( ub(:,:,jk) + z2dt * ua(:,:,jk) ) * umask(:,:,jk)
            va(:,:,jk) = ( vb(:,:,jk) + z2dt * va(:,:,jk) ) * vmask(:,:,jk)
         END DO
      ELSE                                            ! applied on thickness weighted velocity
         DO jk = 1, jpkm1
            ua(:,:,jk) = (          ub(:,:,jk) * fse3u_b(:,:,jk)      &
               &           + z2dt * ua(:,:,jk) * fse3u_n(:,:,jk)  )   &
               &         / fse3u_a(:,:,jk) * umask(:,:,jk)
            va(:,:,jk) = (          vb(:,:,jk) * fse3v_b(:,:,jk)      &
               &           + z2dt * va(:,:,jk) * fse3v_n(:,:,jk)  )   &
               &         / fse3v_a(:,:,jk) * vmask(:,:,jk)
         END DO
      ENDIF
# endif

# if defined key_dynspg_ts
!!gm IF ( lk_dynspg_ts ) THEN ....
      ! Ensure below that barotropic velocities match time splitting estimate
      ! Compute actual transport and replace it with ts estimate at "after" time step
      zue(:,:) = fse3u_a(:,:,1) * ua(:,:,1) * umask(:,:,1)
      zve(:,:) = fse3v_a(:,:,1) * va(:,:,1) * vmask(:,:,1)
      DO jk = 2, jpkm1
         zue(:,:) = zue(:,:) + fse3u_a(:,:,jk) * ua(:,:,jk) * umask(:,:,jk)
         zve(:,:) = zve(:,:) + fse3v_a(:,:,jk) * va(:,:,jk) * vmask(:,:,jk)
      END DO
      DO jk = 1, jpkm1
         ua(:,:,jk) = ( ua(:,:,jk) - zue(:,:) * hur_a(:,:) + ua_b(:,:) ) * umask(:,:,jk)
         va(:,:,jk) = ( va(:,:,jk) - zve(:,:) * hvr_a(:,:) + va_b(:,:) ) * vmask(:,:,jk)
      END DO

      IF (lk_dynspg_ts.AND.(.NOT.ln_bt_fw)) THEN
         ! Remove advective velocity from "now velocities" 
         ! prior to asselin filtering     
         ! In the forward case, this is done below after asselin filtering   
         ! so that asselin contribution is removed at the same time 
         DO jk = 1, jpkm1
            un(:,:,jk) = ( un(:,:,jk) - un_adv(:,:) + un_b(:,:) )*umask(:,:,jk)
            vn(:,:,jk) = ( vn(:,:,jk) - vn_adv(:,:) + vn_b(:,:) )*vmask(:,:,jk)
         END DO  
      ENDIF
!!gm ENDIF
# endif

      ! Update after velocity on domain lateral boundaries
      ! --------------------------------------------------      
      CALL lbc_lnk( ua, 'U', -1. )     !* local domain boundaries
      CALL lbc_lnk( va, 'V', -1. ) 
      !
# if defined key_bdy
      !                                !* BDY open boundaries
      IF( lk_bdy .AND. lk_dynspg_exp ) CALL bdy_dyn( kt )
      IF( lk_bdy .AND. lk_dynspg_ts  ) CALL bdy_dyn( kt, dyn3d_only=.true. )

!!$   Do we need a call to bdy_vol here??
      !
# endif
      !
# if defined key_agrif
      CALL Agrif_dyn( kt )             !* AGRIF zoom boundaries
# endif
#endif

      IF( l_trddyn ) THEN             ! prepare the atf trend computation + some diagnostics
         z1_2dt = 1._wp / (2. * rdt)        ! Euler or leap-frog time step 
         IF( neuler == 0 .AND. kt == nit000 )   z1_2dt = 1._wp / rdt
         !
         !                                  ! Kinetic energy and Conversion
         IF( ln_KE_trd  )   CALL trd_dyn( ua, va, jpdyn_ken, kt )
         !
         IF( ln_dyn_trd ) THEN              ! 3D output: total momentum trends
            zua(:,:,:) = ( ua(:,:,:) - ub(:,:,:) ) * z1_2dt
            zva(:,:,:) = ( va(:,:,:) - vb(:,:,:) ) * z1_2dt
            CALL iom_put( "utrd_tot", zua )        ! total momentum trends, except the asselin time filter
            CALL iom_put( "vtrd_tot", zva )
         ENDIF
         !
         zua(:,:,:) = un(:,:,:)             ! save the now velocity before the asselin filter
         zva(:,:,:) = vn(:,:,:)             ! (caution: there will be a shift by 1 timestep in the
         !                                  !  computation of the asselin filter trends)
      ENDIF

      IF( lk_diavecq ) THEN ! store before and after velocity for thermal wind imbalance trend term

         ua_vecq(:,:,:) = ua(:,:,:)
         va_vecq(:,:,:) = va(:,:,:)

         ub_vecq(:,:,:) = ub(:,:,:)
         vb_vecq(:,:,:) = vb(:,:,:)

      ENDIF

      ! Time filter and swap of dynamics arrays
      ! ------------------------------------------
      IF( neuler == 0 .AND. kt == nit000 ) THEN        !* Euler at first time-step: only swap
         DO jk = 1, jpkm1
            un(:,:,jk) = ua(:,:,jk)                          ! un <-- ua
            vn(:,:,jk) = va(:,:,jk)
         END DO
         IF (lk_vvl) THEN
            DO jk = 1, jpkm1
               fse3t_b(:,:,jk) = fse3t_n(:,:,jk)
               fse3u_b(:,:,jk) = fse3u_n(:,:,jk)
               fse3v_b(:,:,jk) = fse3v_n(:,:,jk)
            ENDDO
         ENDIF
      ELSE                                             !* Leap-Frog : Asselin filter and swap
         !                                ! =============!
         IF( .NOT. lk_vvl ) THEN          ! Fixed volume !
            !                             ! =============!
            DO jk = 1, jpkm1                              
               DO jj = 1, jpj
                  DO ji = 1, jpi    
                     zuf = un(ji,jj,jk) + atfp * ( ub(ji,jj,jk) - 2._wp * un(ji,jj,jk) + ua(ji,jj,jk) )
                     zvf = vn(ji,jj,jk) + atfp * ( vb(ji,jj,jk) - 2._wp * vn(ji,jj,jk) + va(ji,jj,jk) )
                     !
                     ub(ji,jj,jk) = zuf                      ! ub <-- filtered velocity
                     vb(ji,jj,jk) = zvf
                     un(ji,jj,jk) = ua(ji,jj,jk)             ! un <-- ua
                     vn(ji,jj,jk) = va(ji,jj,jk)
                  END DO
               END DO
            END DO
            !                             ! ================!
         ELSE                             ! Variable volume !
            !                             ! ================!
            ! Before scale factor at t-points
            ! (used as a now filtered scale factor until the swap)
            ! ----------------------------------------------------
            IF (lk_dynspg_ts.AND.ln_bt_fw) THEN
               ! No asselin filtering on thicknesses if forward time splitting
                  fse3t_b(:,:,:) = fse3t_n(:,:,:)
            ELSE
               fse3t_b(:,:,:) = fse3t_n(:,:,:) + atfp * ( fse3t_b(:,:,:) - 2._wp * fse3t_n(:,:,:) + fse3t_a(:,:,:) )
               ! Add volume filter correction: compatibility with tracer advection scheme
               ! => time filter + conservation correction (only at the first level)
               IF ( nn_isf == 0) THEN   ! if no ice shelf melting
                  fse3t_b(:,:,1) = fse3t_b(:,:,1) - atfp * rdt * r1_rau0 * ( emp_b(:,:) - emp(:,:) &
                                 &                                          -rnf_b(:,:) + rnf(:,:) ) * tmask(:,:,1)
               ELSE                     ! if ice shelf melting
                  DO jj = 1,jpj
                     DO ji = 1,jpi
                        jk = mikt(ji,jj)
                        fse3t_b(ji,jj,jk) = fse3t_b(ji,jj,jk) - atfp * rdt * r1_rau0                       &
                                          &                          * ( (emp_b(ji,jj)    - emp(ji,jj)   ) &
                                          &                            - (rnf_b(ji,jj)    - rnf(ji,jj)   ) &
                                          &                            + (fwfisf_b(ji,jj) - fwfisf(ji,jj)) ) * tmask(ji,jj,jk)
                     END DO
                  END DO
               END IF
            ENDIF
            !
            IF( ln_dynadv_vec ) THEN
               ! Before scale factor at (u/v)-points
               ! -----------------------------------
               CALL dom_vvl_interpol( fse3t_b(:,:,:), fse3u_b(:,:,:), 'U' )
               CALL dom_vvl_interpol( fse3t_b(:,:,:), fse3v_b(:,:,:), 'V' )
               ! Leap-Frog - Asselin filter and swap: applied on velocity
               ! -----------------------------------
               DO jk = 1, jpkm1
                  DO jj = 1, jpj
                     DO ji = 1, jpi
                        zuf = un(ji,jj,jk) + atfp * ( ub(ji,jj,jk) - 2._wp * un(ji,jj,jk) + ua(ji,jj,jk) )
                        zvf = vn(ji,jj,jk) + atfp * ( vb(ji,jj,jk) - 2._wp * vn(ji,jj,jk) + va(ji,jj,jk) )
                        !
                        ub(ji,jj,jk) = zuf                      ! ub <-- filtered velocity
                        vb(ji,jj,jk) = zvf
                        un(ji,jj,jk) = ua(ji,jj,jk)             ! un <-- ua
                        vn(ji,jj,jk) = va(ji,jj,jk)
                     END DO
                  END DO
               END DO
               !
            ELSE
               ! Temporary filtered scale factor at (u/v)-points (will become before scale factor)
               !------------------------------------------------
               CALL dom_vvl_interpol( fse3t_b(:,:,:), ze3u_f, 'U' )
               CALL dom_vvl_interpol( fse3t_b(:,:,:), ze3v_f, 'V' )
               ! Leap-Frog - Asselin filter and swap: applied on thickness weighted velocity
               ! -----------------------------------             ===========================
               DO jk = 1, jpkm1
                  DO jj = 1, jpj
                     DO ji = 1, jpi                  
                        zue3a = ua(ji,jj,jk) * fse3u_a(ji,jj,jk)
                        zve3a = va(ji,jj,jk) * fse3v_a(ji,jj,jk)
                        zue3n = un(ji,jj,jk) * fse3u_n(ji,jj,jk)
                        zve3n = vn(ji,jj,jk) * fse3v_n(ji,jj,jk)
                        zue3b = ub(ji,jj,jk) * fse3u_b(ji,jj,jk)
                        zve3b = vb(ji,jj,jk) * fse3v_b(ji,jj,jk)
                        !
                        zuf = ( zue3n + atfp * ( zue3b - 2._wp * zue3n  + zue3a ) ) / ze3u_f(ji,jj,jk)
                        zvf = ( zve3n + atfp * ( zve3b - 2._wp * zve3n  + zve3a ) ) / ze3v_f(ji,jj,jk)
                        !
                        ub(ji,jj,jk) = zuf                     ! ub <-- filtered velocity
                        vb(ji,jj,jk) = zvf
                        un(ji,jj,jk) = ua(ji,jj,jk)            ! un <-- ua
                        vn(ji,jj,jk) = va(ji,jj,jk)
                     END DO
                  END DO
               END DO
               fse3u_b(:,:,1:jpkm1) = ze3u_f(:,:,1:jpkm1)      ! e3u_b <-- filtered scale factor
               fse3v_b(:,:,1:jpkm1) = ze3v_f(:,:,1:jpkm1)
            ENDIF
            !
         ENDIF
         !
         IF (lk_dynspg_ts.AND.ln_bt_fw) THEN
            ! Revert "before" velocities to time split estimate
            ! Doing it here also means that asselin filter contribution is removed  
            zue(:,:) = fse3u_b(:,:,1) * ub(:,:,1) * umask(:,:,1)
            zve(:,:) = fse3v_b(:,:,1) * vb(:,:,1) * vmask(:,:,1)    
            DO jk = 2, jpkm1
               zue(:,:) = zue(:,:) + fse3u_b(:,:,jk) * ub(:,:,jk) * umask(:,:,jk)
               zve(:,:) = zve(:,:) + fse3v_b(:,:,jk) * vb(:,:,jk) * vmask(:,:,jk)    
            END DO
            DO jk = 1, jpkm1
               ub(:,:,jk) = ub(:,:,jk) - (zue(:,:) * hur(:,:) - un_b(:,:)) * umask(:,:,jk)
               vb(:,:,jk) = vb(:,:,jk) - (zve(:,:) * hvr(:,:) - vn_b(:,:)) * vmask(:,:,jk)
            END DO
         ENDIF
         !
      ENDIF ! neuler =/0
      !
      ! Set "now" and "before" barotropic velocities for next time step:
      ! JC: Would be more clever to swap variables than to make a full vertical
      ! integration
      !
      !
      IF (lk_vvl) THEN
         hu_b(:,:) = 0.
         hv_b(:,:) = 0.
         DO jk = 1, jpkm1
            hu_b(:,:) = hu_b(:,:) + fse3u_b(:,:,jk) * umask(:,:,jk)
            hv_b(:,:) = hv_b(:,:) + fse3v_b(:,:,jk) * vmask(:,:,jk)
         END DO
         hur_b(:,:) = umask_i(:,:) / ( hu_b(:,:) + 1._wp - umask_i(:,:) )
         hvr_b(:,:) = vmask_i(:,:) / ( hv_b(:,:) + 1._wp - vmask_i(:,:) )
      ENDIF
      !
      un_b(:,:) = 0._wp ; vn_b(:,:) = 0._wp
      ub_b(:,:) = 0._wp ; vb_b(:,:) = 0._wp
      !
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               un_b(ji,jj) = un_b(ji,jj) + fse3u_a(ji,jj,jk) * un(ji,jj,jk) * umask(ji,jj,jk)
               vn_b(ji,jj) = vn_b(ji,jj) + fse3v_a(ji,jj,jk) * vn(ji,jj,jk) * vmask(ji,jj,jk)
               !
               ub_b(ji,jj) = ub_b(ji,jj) + fse3u_b(ji,jj,jk) * ub(ji,jj,jk) * umask(ji,jj,jk)
               vb_b(ji,jj) = vb_b(ji,jj) + fse3v_b(ji,jj,jk) * vb(ji,jj,jk) * vmask(ji,jj,jk)
            END DO
         END DO
      END DO
      !
      !
      un_b(:,:) = un_b(:,:) * hur_a(:,:)
      vn_b(:,:) = vn_b(:,:) * hvr_a(:,:)
      ub_b(:,:) = ub_b(:,:) * hur_b(:,:)
      vb_b(:,:) = vb_b(:,:) * hvr_b(:,:)
      !
      !

      IF( l_trddyn ) THEN                ! 3D output: asselin filter trends on momentum
         zua(:,:,:) = ( ub(:,:,:) - zua(:,:,:) ) * z1_2dt
         zva(:,:,:) = ( vb(:,:,:) - zva(:,:,:) ) * z1_2dt
         CALL trd_dyn( zua, zva, jpdyn_atf, kt )
      ENDIF
      !
      IF(ln_ctl)   CALL prt_ctl( tab3d_1=un, clinfo1=' nxt  - Un: ', mask1=umask,   &
         &                       tab3d_2=vn, clinfo2=' Vn: '       , mask2=vmask )
      ! 
      CALL wrk_dealloc( jpi,jpj,jpk,  ze3u_f, ze3v_f, zua, zva )
      IF( lk_dynspg_ts )   CALL wrk_dealloc( jpi,jpj, zue, zve )
      !
      IF( nn_timing == 1 )  CALL timing_stop('dyn_nxt')
      !
   END SUBROUTINE dyn_nxt

   !!=========================================================================
END MODULE dynnxt
