MODULE step
   !!======================================================================
   !!                       ***  MODULE step  ***
   !! Time-stepping    : manager of the ocean, tracer and ice time stepping
   !!======================================================================
   !! History :  OPA  !  1991-03  (G. Madec)  Original code
   !!             -   !  1991-11  (G. Madec)
   !!             -   !  1992-06  (M. Imbard)  add a first output record
   !!             -   !  1996-04  (G. Madec)  introduction of dynspg
   !!             -   !  1996-04  (M.A. Foujols)  introduction of passive tracer
   !!            8.0  !  1997-06  (G. Madec)  new architecture of call
   !!            8.2  !  1997-06  (G. Madec, M. Imbard, G. Roullet)  free surface
   !!             -   !  1999-02  (G. Madec, N. Grima)  hpg implicit
   !!             -   !  2000-07  (J-M Molines, M. Imbard)  Open Bondary Conditions
   !!   NEMO     1.0  !  2002-06  (G. Madec)  free form, suppress macro-tasking
   !!             -   !  2004-08  (C. Talandier) New trends organization
   !!             -   !  2005-01  (C. Ethe) Add the KPP closure scheme
   !!             -   !  2005-11  (G. Madec)  Reorganisation of tra and dyn calls
   !!             -   !  2006-01  (L. Debreu, C. Mazauric)  Agrif implementation
   !!             -   !  2006-07  (S. Masson)  restart using iom
   !!            3.2  !  2009-02  (G. Madec, R. Benshila)  reintroduicing z*-coordinate
   !!             -   !  2009-06  (S. Masson, G. Madec)  TKE restart compatible with key_cpl
   !!            3.3  !  2010-05  (K. Mogensen, A. Weaver, M. Martin, D. Lea) Assimilation interface
   !!             -   !  2010-10  (C. Ethe, G. Madec) reorganisation of initialisation phase + merge TRC-TRA
   !!            3.4  !  2011-04  (G. Madec, C. Ethe) Merge of dtatem and dtasal
   !!                 !  2012-07  (J. Simeon, G. Madec, C. Ethe) Online coarsening of outputs
   !!            3.7  !  2014-04  (F. Roquet, G. Madec) New equations of state
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   stp             : OPA system time-stepping
   !!----------------------------------------------------------------------
   USE step_oce         ! time stepping definition modules
   USE iom
   USE diavecq
   IMPLICIT NONE
   PRIVATE

   PUBLIC   stp   ! called by opa.F90

   !! * Substitutions
#  include "domzgr_substitute.h90"
!!gm   #  include "zdfddm_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.7 , NEMO Consortium (2014)
   !! $Id: step.F90 6204 2016-01-04 13:47:06Z cetlod $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

#if defined key_agrif
   RECURSIVE SUBROUTINE stp( )
      INTEGER             ::   kstp   ! ocean time-step index
#else
   SUBROUTINE stp( kstp )
      INTEGER, INTENT(in) ::   kstp   ! ocean time-step index
#endif
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE stp  ***
      !!
      !! ** Purpose : - Time stepping of OPA (momentum and active tracer eqs.)
      !!              - Time stepping of LIM (dynamic and thermodynamic eqs.)
      !!              - Tme stepping  of TRC (passive tracer eqs.)
      !!
      !! ** Method  : -1- Update forcings and data
      !!              -2- Update ocean physics
      !!              -3- Compute the t and s trends
      !!              -4- Update t and s
      !!              -5- Compute the momentum trends
      !!              -6- Update the horizontal velocity
      !!              -7- Compute the diagnostics variables (rd,N2, div,cur,w)
      !!              -8- Outputs and diagnostics
      !!----------------------------------------------------------------------
      INTEGER ::   jk       ! dummy loop indice
      INTEGER ::   indic    ! error indicator if < 0
      INTEGER ::   kcall    ! optional integer argument (dom_vvl_sf_nxt)
      REAL(wp)  ::   zstp
      !! ---------------------------------------------------------------------

#if defined key_agrif
      kstp = nit000 + Agrif_Nb_Step()
      IF ( lk_agrif_debug ) THEN
         IF ( Agrif_Root() .and. lwp) Write(*,*) '---'
         IF (lwp) Write(*,*) 'Grid Number',Agrif_Fixed(),' time step ',kstp, 'int tstep',Agrif_NbStepint()
      ENDIF

      IF ( kstp == (nit000 + 1) ) lk_agrif_fstep = .FALSE.


# if defined key_iomput
      IF( Agrif_Nbstepint() == 0 )   CALL iom_swap( cxios_context )
# endif
#endif

#if defined key_diavecq   
      zstp=kstp-nit000+1
      IF ( MOD(zstp,86400._wp/rdt) == 0.)  THEN
        l_diavecq_out = .TRUE.   
      ELSE
       l_diavecq_out = .FALSE.   
!       l_diavecq_out = .TRUE.   
      ENDIF
#endif

                             indic = 0           ! reset to no error condition
      IF( kstp == nit000 ) THEN
         ! must be done after nemo_init for AGRIF+XIOS+OASIS
                      CALL iom_init(      cxios_context          )  ! iom_put initialization
         IF( ln_crs ) CALL iom_init( TRIM(cxios_context)//"_crs" )  ! initialize context for coarse grid
      ENDIF

      IF( kstp /= nit000 )   CALL day( kstp )         ! Calendar (day was already called at nit000 in day_init)
                             CALL iom_setkt( kstp - nit000 + 1,      cxios_context          )   ! tell iom we are at time step kstp
      IF( ln_crs     )       CALL iom_setkt( kstp - nit000 + 1, TRIM(cxios_context)//"_crs" )   ! tell iom we are at time step kstp

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Update data, open boundaries, surface boundary condition (including sea-ice)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      IF( lk_tide    )   CALL sbc_tide( kstp )
                         CALL sbc    ( kstp )         ! Sea Boundary Condition (including sea-ice)
      IF( lk_bdy     )  THEN
         IF( ln_apr_dyn) CALL sbc_apr( kstp )   ! bdy_dta needs ssh_ib 
                         CALL bdy_dta ( kstp, time_offset=+1 )   ! update dynamic & tracer data at open boundaries
      ENDIF
!                         CALL sbc    ( kstp )         ! Sea Boundary Condition (including sea-ice)
                                                      ! clem: moved here for bdy ice purpose
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Update stochastic parameters and random T/S fluctuations
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
       IF( ln_sto_eos ) CALL sto_par( kstp )          ! Stochastic parameters
       IF( ln_sto_eos ) CALL sto_pts( tsn  )          ! Random T/S fluctuations

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Ocean physics update                (ua, va, tsa used as workspace)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !  THERMODYNAMICS
                         CALL eos_rab( tsb, rab_b )       ! before local thermal/haline expension ratio at T-points
                         CALL eos_rab( tsn, rab_n )       ! now    local thermal/haline expension ratio at T-points
                         CALL bn2    ( tsb, rab_b, rn2b ) ! before Brunt-Vaisala frequency
                         CALL bn2    ( tsn, rab_n, rn2  ) ! now    Brunt-Vaisala frequency
      !
      !  VERTICAL PHYSICS
                         CALL zdf_bfr( kstp )         ! bottom friction (if quadratic)
      !                                               ! Vertical eddy viscosity and diffusivity coefficients
      IF( lk_zdfric  )   CALL zdf_ric( kstp )            ! Richardson number dependent Kz
      IF( lk_zdftke  )   CALL zdf_tke( kstp )            ! TKE closure scheme for Kz
      IF( lk_zdfgls  )   CALL zdf_gls( kstp )            ! GLS closure scheme for Kz
      IF( lk_zdfkpp  )   CALL zdf_kpp( kstp )            ! KPP closure scheme for Kz
      IF( lk_zdfcst  ) THEN                              ! Constant Kz (reset avt, avm[uv] to the background value)
         avt (:,:,:) = rn_avt0 * wmask (:,:,:)
         avmu(:,:,:) = rn_avm0 * wumask(:,:,:)
         avmv(:,:,:) = rn_avm0 * wvmask(:,:,:)
      ENDIF
      IF( ln_rnf_mouth ) THEN                         ! increase diffusivity at rivers mouths
         DO jk = 2, nkrnf   ;   avt(:,:,jk) = avt(:,:,jk) + 2.e0 * rn_avt_rnf * rnfmsk(:,:) * tmask(:,:,jk)   ;   END DO
      ENDIF
      IF( ln_zdfevd  )   CALL zdf_evd( kstp )         ! enhanced vertical eddy diffusivity

      IF( lk_zdftmx  )   CALL zdf_tmx( kstp )         ! tidal vertical mixing

      IF( lk_zdfddm .AND. .NOT. lk_zdfkpp )   &
         &               CALL zdf_ddm( kstp )         ! double diffusive mixing

                         CALL zdf_mxl( kstp )         ! mixed layer depth

                                                      ! write TKE or GLS information in the restart file
      IF( lrst_oce .AND. lk_zdftke )   CALL tke_rst( kstp, 'WRITE' )
      IF( lrst_oce .AND. lk_zdfgls )   CALL gls_rst( kstp, 'WRITE' )
      !
      !  LATERAL  PHYSICS
      !
      IF( lk_ldfslp ) THEN                            ! slope of lateral mixing
                         CALL eos( tsb, rhd, gdept_0(:,:,:) )               ! before in situ density
         IF( ln_zps .AND. .NOT. ln_isfcav)                               &
            &            CALL zps_hde    ( kstp, jpts, tsb, gtsu, gtsv,  &  ! Partial steps: before horizontal gradient
            &                                          rhd, gru , grv    )  ! of t, s, rd at the last ocean level
         IF( ln_zps .AND.       ln_isfcav)                               &
            &            CALL zps_hde_isf( kstp, jpts, tsb, gtsu, gtsv,  &    ! Partial steps for top cell (ISF)
            &                                          rhd, gru , grv , aru , arv , gzu , gzv , ge3ru , ge3rv ,   &
            &                                   gtui, gtvi, grui, grvi, arui, arvi, gzui, gzvi, ge3rui, ge3rvi    ) ! of t, s, rd at the first ocean level
         IF( ln_traldf_grif ) THEN                           ! before slope for Griffies operator
                         CALL ldf_slp_grif( kstp )
         ELSE
                         CALL ldf_slp( kstp, rhd, rn2b )     ! before slope for Madec operator
         ENDIF
      ENDIF
#if defined key_traldf_c2d
      IF( lk_traldf_eiv )   CALL ldf_eiv( kstp )      ! eddy induced velocity coefficient
#endif
#if defined key_traldf_c3d && defined key_traldf_smag
                          CALL ldf_tra_smag( kstp )      ! eddy induced velocity coefficient
#  endif
#if defined key_dynldf_c3d && defined key_dynldf_smag
                          CALL ldf_dyn_smag( kstp )      ! eddy induced velocity coefficient
#  endif

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !  Ocean dynamics : hdiv, rot, ssh, e3, wn
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                         CALL ssh_nxt       ( kstp )  ! after ssh (includes call to div_cur)
      IF( lk_vvl     )   CALL dom_vvl_sf_nxt( kstp )  ! after vertical scale factors 
                         CALL wzv           ( kstp )  ! now cross-level velocity 

      IF( lk_dynspg_ts ) THEN 
          ! In case the time splitting case, update almost all momentum trends here:
          ! Note that the computation of vertical velocity above, hence "after" sea level
          ! is necessary to compute momentum advection for the rhs of barotropic loop:
                            CALL eos    ( tsn, rhd, rhop, fsdept_n(:,:,:) ) ! now in situ density for hpg computation
            IF( ln_zps .AND. .NOT. ln_isfcav)                               &
               &            CALL zps_hde    ( kstp, jpts, tsn, gtsu, gtsv,  &    ! Partial steps: before horizontal gradient
               &                                          rhd, gru , grv    )  ! of t, s, rd at the last ocean level
            IF( ln_zps .AND.       ln_isfcav)                               &
               &            CALL zps_hde_isf( kstp, jpts, tsn, gtsu, gtsv,  &    ! Partial steps for top cell (ISF)
               &                                          rhd, gru , grv , aru , arv , gzu , gzv , ge3ru , ge3rv ,   &
               &                                   gtui, gtvi, grui, grvi, arui, arvi, gzui, gzvi, ge3rui, ge3rvi    ) ! of t, s, rd at the last ocean level

                                  ua(:,:,:) = 0.e0             ! set dynamics trends to zero
                                  va(:,:,:) = 0.e0
          IF(  lk_asminc .AND. ln_asmiau .AND. &
             & ln_dyninc       )  CALL dyn_asm_inc  ( kstp )   ! apply dynamics assimilation increment
          IF( ln_neptsimp )       CALL dyn_nept_cor ( kstp )   ! subtract Neptune velocities (simplified)
          IF( lk_bdy           )  CALL bdy_dyn3d_dmp( kstp )   ! bdy damping trends
                                  CALL dyn_adv      ( kstp )   ! advection (vector or flux form)
                                  CALL dyn_vor      ( kstp )   ! vorticity term including Coriolis
                                  CALL dyn_ldf      ( kstp )   ! lateral mixing
          IF( ln_neptsimp )       CALL dyn_nept_cor ( kstp )   ! add Neptune velocities (simplified)
#if defined key_agrif
          IF(.NOT. Agrif_Root())  CALL Agrif_Sponge_dyn        ! momentum sponge
#endif
                                  CALL dyn_hpg( kstp )         ! horizontal gradient of Hydrostatic pressure
                                  CALL dyn_spg( kstp, indic )  ! surface pressure gradient

                                  ua_sv(:,:,:) = ua(:,:,:)     ! Save trends (barotropic trend has been fully updated at this stage)
                                  va_sv(:,:,:) = va(:,:,:)

                                  CALL div_cur( kstp )         ! Horizontal divergence & Relative vorticity (2nd call in time-split case)
          IF( lk_vvl     )        CALL dom_vvl_sf_nxt( kstp, kcall=2 )  ! after vertical scale factors (update depth average component)
                                  CALL wzv           ( kstp )  ! now cross-level velocity 
      ENDIF

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! diagnostics and outputs             (ua, va, tsa used as workspace)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      IF( lk_floats  )      CALL flo_stp( kstp )         ! drifting Floats
      IF( lk_diahth  )      CALL dia_hth( kstp )         ! Thermocline depth (20 degres isotherm depth)
      IF( .NOT. ln_cpl )    CALL dia_fwb( kstp )         ! Fresh water budget diagnostics
      IF( lk_diadct  )      CALL dia_dct( kstp )         ! Transports
      IF( lk_diaar5  )      CALL dia_ar5( kstp )         ! ar5 diag
      IF( lk_diavecq )      CALL dia_vecq( kstp )        ! vector Q components diags
      IF( lk_diaharm )      CALL dia_harm( kstp )        ! Tidal harmonic analysis
                            CALL dia_wri( kstp )         ! ocean model: outputs
      !
      IF( ln_crs     )      CALL crs_fld( kstp )         ! ocean model: online field coarsening & output

#if defined key_top
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Passive Tracer Model
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                         CALL trc_stp( kstp )         ! time-stepping
#endif


      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Active tracers                              (ua, va used as workspace)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                             tsa(:,:,:,:) = 0.e0            ! set tracer trends to zero

      IF(  lk_asminc .AND. ln_asmiau .AND. &
         & ln_trainc     )   CALL tra_asm_inc( kstp )       ! apply tracer assimilation increment
                             CALL tra_sbc    ( kstp )       ! surface boundary condition
      IF( ln_traqsr      )   CALL tra_qsr    ( kstp )       ! penetrative solar radiation qsr
      IF( ln_trabbc      )   CALL tra_bbc    ( kstp )       ! bottom heat flux
      IF( lk_trabbl      )   CALL tra_bbl    ( kstp )       ! advective (and/or diffusive) bottom boundary layer scheme
      IF( ln_tradmp      )   CALL tra_dmp    ( kstp )       ! internal damping trends
      IF( lk_bdy         )   CALL bdy_tra_dmp( kstp )       ! bdy damping trends
                             CALL tra_adv    ( kstp )       ! horizontal & vertical advection
      IF( lk_zdfkpp      )   CALL tra_kpp    ( kstp )       ! KPP non-local tracer fluxes
                             CALL tra_ldf    ( kstp )       ! lateral mixing

      IF( ln_diaptr      )   CALL dia_ptr                   ! Poleward adv/ldf TRansports diagnostics

#if defined key_agrif
      IF(.NOT. Agrif_Root()) CALL Agrif_Sponge_tra          ! tracers sponge
#endif
                             CALL tra_zdf    ( kstp )       ! vertical mixing and after tracer fields

      IF( ln_dynhpg_imp  ) THEN                             ! semi-implicit hpg (time stepping then eos)
         IF( ln_zdfnpc   )   CALL tra_npc( kstp )                ! update after fields by non-penetrative convection
                             CALL tra_nxt( kstp )                ! tracer fields at next time step
                             CALL eos    ( tsa, rhd, rhop, fsdept_n(:,:,:) )  ! Time-filtered in situ density for hpg computation
            IF( ln_zps .AND. .NOT. ln_isfcav)                                &
               &             CALL zps_hde    ( kstp, jpts, tsa, gtsu, gtsv,  &    ! Partial steps: before horizontal gradient
               &                                           rhd, gru , grv    )  ! of t, s, rd at the last ocean level
            IF( ln_zps .AND.       ln_isfcav)                                &
               &             CALL zps_hde_isf( kstp, jpts, tsa, gtsu, gtsv,  &    ! Partial steps for top cell (ISF)
               &                                           rhd, gru , grv , aru , arv , gzu , gzv , ge3ru , ge3rv ,   &
               &                                    gtui, gtvi, grui, grvi, arui, arvi, gzui, gzvi, ge3rui, ge3rvi    ) ! of t, s, rd at the last ocean level
      ELSE                                                  ! centered hpg  (eos then time stepping)
         IF ( .NOT. lk_dynspg_ts ) THEN                     ! eos already called in time-split case
                             CALL eos    ( tsn, rhd, rhop, fsdept_n(:,:,:) )  ! now in situ density for hpg computation
         IF( ln_zps .AND. .NOT. ln_isfcav)                                   &
               &             CALL zps_hde    ( kstp, jpts, tsn, gtsu, gtsv,  &    ! Partial steps: before horizontal gradient
               &                                           rhd, gru , grv    )  ! of t, s, rd at the last ocean level
         IF( ln_zps .AND.       ln_isfcav)                                   & 
               &             CALL zps_hde_isf( kstp, jpts, tsn, gtsu, gtsv,  &    ! Partial steps for top cell (ISF)
               &                                           rhd, gru , grv , aru , arv , gzu , gzv , ge3ru , ge3rv ,   &
               &                                    gtui, gtvi, grui, grvi, arui, arvi, gzui, gzvi, ge3rui, ge3rvi    ) ! of t, s, rd at the last ocean level
         ENDIF
         IF( ln_zdfnpc   )   CALL tra_npc( kstp )                ! update after fields by non-penetrative convection
                             CALL tra_nxt( kstp )                ! tracer fields at next time step
      ENDIF

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Dynamics                                    (tsa used as workspace)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      IF( lk_dynspg_ts   )  THEN
                                                             ! revert to previously computed momentum tendencies
                                                             ! (not using ua, va as temporary arrays during tracers' update could avoid that)
                               ua(:,:,:) = ua_sv(:,:,:)
                               va(:,:,:) = va_sv(:,:,:)
                                                             ! Revert now divergence and rotational to previously computed ones 
                                                             !(needed because of the time swap in div_cur, at the beginning of each time step)
                               hdivn(:,:,:) = hdivb(:,:,:)
                               rotn(:,:,:)  = rotb(:,:,:) 

                               CALL dyn_bfr( kstp )         ! bottom friction
                               CALL dyn_zdf( kstp )         ! vertical diffusion
      ELSE
                               ua(:,:,:) = 0.e0             ! set dynamics trends to zero
                               va(:,:,:) = 0.e0

        IF(  lk_asminc .AND. ln_asmiau .AND. &
           & ln_dyninc      )  CALL dyn_asm_inc( kstp )     ! apply dynamics assimilation increment
        IF( ln_bkgwri )        CALL asm_bkg_wri( kstp )     ! output background fields
        IF( ln_neptsimp )      CALL dyn_nept_cor( kstp )    ! subtract Neptune velocities (simplified)
        IF( lk_bdy          )  CALL bdy_dyn3d_dmp(kstp )    ! bdy damping trends
                               CALL dyn_adv( kstp )         ! advection (vector or flux form)
                               CALL dyn_vor( kstp )         ! vorticity term including Coriolis
                               CALL dyn_ldf( kstp )         ! lateral mixing
        IF( ln_neptsimp )      CALL dyn_nept_cor( kstp )    ! add Neptune velocities (simplified)
#if defined key_agrif
        IF(.NOT. Agrif_Root()) CALL Agrif_Sponge_dyn        ! momemtum sponge
#endif
                               CALL dyn_hpg( kstp )         ! horizontal gradient of Hydrostatic pressure
                               CALL dyn_bfr( kstp )         ! bottom friction
                               CALL dyn_zdf( kstp )         ! vertical diffusion
                               CALL dyn_spg( kstp, indic )  ! surface pressure gradient
      ENDIF
                               CALL dyn_nxt( kstp )         ! lateral velocity at next time step

                               CALL ssh_swp( kstp )         ! swap of sea surface height
      IF( lk_vvl           )   CALL dom_vvl_sf_swp( kstp )  ! swap of vertical scale factors
      !
      IF( lrst_oce         )   CALL rst_write( kstp )       ! write output ocean restart file

#if defined key_agrif
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! AGRIF
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<      
                               CALL Agrif_Integrate_ChildGrids( stp )  

      IF ( Agrif_NbStepint().EQ.0 ) THEN
                               CALL Agrif_Update_Tra()      ! Update active tracers
                               CALL Agrif_Update_Dyn()      ! Update momentum
      ENDIF
#endif
      IF( ln_diahsb        )   CALL dia_hsb( kstp )         ! - ML - global conservation diagnostics
      IF( lk_diaobs  )         CALL dia_obs( kstp )         ! obs-minus-model (assimilation) diagnostics (call after dynamics update)

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Control
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                               CALL stp_ctl( kstp, indic )
      IF( indic < 0        )   THEN
                               CALL ctl_stop( 'step: indic < 0' )
                               CALL dia_wri_state( 'output.abort', kstp )
      ENDIF
      IF( kstp == nit000   )   THEN
                 CALL iom_close( numror )     ! close input  ocean restart file
         IF(lwm) CALL FLUSH    ( numond )     ! flush output namelist oce
         IF( lwm.AND.numoni /= -1 ) CALL FLUSH    ( numoni )     ! flush output namelist ice
      ENDIF

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Coupled mode
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      IF( lk_oasis         )   CALL sbc_cpl_snd( kstp )     ! coupled mode : field exchanges
      !
#if defined key_iomput
      IF( kstp == nitend .OR. indic < 0 ) THEN 
                      CALL iom_context_finalize(      cxios_context          ) ! needed for XIOS+AGRIF
         IF( ln_crs ) CALL iom_context_finalize( trim(cxios_context)//"_crs" ) ! 
      ENDIF
#endif
      !
      IF( nn_timing == 1 .AND.  kstp == nit000  )   CALL timing_reset
      !     
      !
   END SUBROUTINE stp

   !!======================================================================
END MODULE step
