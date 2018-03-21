MODULE tranxt
   !!======================================================================
   !!                       ***  MODULE  tranxt  ***
   !! Ocean active tracers:  time stepping on temperature and salinity
   !!======================================================================
   !! History :  OPA  !  1991-11  (G. Madec)  Original code
   !!            7.0  !  1993-03  (M. Guyon)  symetrical conditions
   !!            8.0  !  1996-02  (G. Madec & M. Imbard)  opa release 8.0
   !!             -   !  1996-04  (A. Weaver)  Euler forward step
   !!            8.2  !  1999-02  (G. Madec, N. Grima)  semi-implicit pressure grad.
   !!  NEMO      1.0  !  2002-08  (G. Madec)  F90: Free form and module
   !!             -   !  2002-11  (C. Talandier, A-M Treguier) Open boundaries
   !!             -   !  2005-04  (C. Deltel) Add Asselin trend in the ML budget
   !!            2.0  !  2006-02  (L. Debreu, C. Mazauric) Agrif implementation
   !!            3.0  !  2008-06  (G. Madec)  time stepping always done in trazdf
   !!            3.1  !  2009-02  (G. Madec, R. Benshila)  re-introduce the vvl option
   !!            3.3  !  2010-04  (M. Leclair, G. Madec)  semi-implicit hpg with asselin filter + modified LF-RA
   !!             -   !  2010-05  (C. Ethe, G. Madec)  merge TRC-TRA
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   tra_nxt       : time stepping on tracers
   !!   tra_nxt_fix   : time stepping on tracers : fixed    volume case
   !!   tra_nxt_vvl   : time stepping on tracers : variable volume case
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers variables
   USE iom
   USE dom_oce         ! ocean space and time domain variables 
   USE sbc_oce         ! surface boundary condition: ocean
   USE sbcrnf          ! river runoffs
   USE sbcisf          ! ice shelf melting/freezing
   USE zdf_oce         ! ocean vertical mixing
   USE domvvl          ! variable volume
   USE dynspg_oce      ! surface     pressure gradient variables
   USE dynhpg          ! hydrostatic pressure gradient 
   USE trd_oce         ! trends: ocean variables
   USE trdtra          ! trends manager: tracers 
   USE diavecq         ! diags for Q vector
   USE traqsr          ! penetrative solar radiation (needed for nksr)
   USE phycst          ! physical constant
   USE ldftra_oce      ! lateral physics on tracers
   USE bdy_oce         ! BDY open boundary condition variables
   USE bdytra          ! open boundary condition (bdy_tra routine)
   !
   USE in_out_manager  ! I/O manager
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE prtctl          ! Print control
   USE wrk_nemo        ! Memory allocation
   USE timing          ! Timing
#if defined key_agrif
   USE agrif_opa_interp
#endif

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tra_nxt       ! routine called by step.F90
   PUBLIC   tra_nxt_fix   ! to be used in trcnxt
   PUBLIC   tra_nxt_vvl   ! to be used in trcnxt

   REAL(wp) ::   rbcp   ! Brown & Campana parameters for semi-implicit hpg

   !! * Substitutions
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO-Consortium (2010) 
   !! $Id: tranxt.F90 6204 2016-01-04 13:47:06Z cetlod $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE tra_nxt( kt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE tranxt  ***
      !!
      !! ** Purpose :   Apply the boundary condition on the after temperature  
      !!             and salinity fields, achieved the time stepping by adding
      !!             the Asselin filter on now fields and swapping the fields.
      !! 
      !! ** Method  :   At this stage of the computation, ta and sa are the 
      !!             after temperature and salinity as the time stepping has
      !!             been performed in trazdf_imp or trazdf_exp module.
      !!
      !!              - Apply lateral boundary conditions on (ta,sa) 
      !!             at the local domain   boundaries through lbc_lnk call, 
      !!             at the one-way open boundaries (lk_bdy=T), 
      !!             at the AGRIF zoom   boundaries (lk_agrif=T)
      !!
      !!              - Update lateral boundary conditions on AGRIF children
      !!             domains (lk_agrif=T)
      !!
      !! ** Action  : - (tb,sb) and (tn,sn) ready for the next time step
      !!              - (ta,sa) time averaged (t,s)   (ln_dynhpg_imp = T)
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt    ! ocean time-step index
      !!
      INTEGER  ::   jk, jn    ! dummy loop indices
      REAL(wp) ::   zfact, z1_2dt, zstp     ! local scalars
      REAL(wp), POINTER, DIMENSION(:,:,:) ::  ztrdt, ztrds
      REAL(wp), POINTER, DIMENSION(:,:,:) ::  zugb, zvgb, zuga, zvga
      REAL(wp), POINTER, DIMENSION(:,:,:) ::  zthub, zthvb, zthua, zthva
      REAL(wp), POINTER, DIMENSION(:,:,:) ::  zdtthu, zdtthv
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start( 'tra_nxt')
      zstp=kt-nit000+1

      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_nxt : achieve the time stepping by Asselin filter and array swap'
         IF(lwp) WRITE(numout,*) '~~~~~~~'
         !
         rbcp = 0.25_wp * (1._wp + atfp) * (1._wp + atfp) * ( 1._wp - atfp)      ! Brown & Campana parameter for semi-implicit hpg
      ENDIF

      ! Update after tracer on domain lateral boundaries
      ! 
#if defined key_agrif
      CALL Agrif_tra                     ! AGRIF zoom boundaries
#endif
      !
      CALL lbc_lnk( tsa(:,:,:,jp_tem), 'T', 1._wp )      ! local domain boundaries  (T-point, unchanged sign)
      CALL lbc_lnk( tsa(:,:,:,jp_sal), 'T', 1._wp )
      !
#if defined key_bdy 
      IF( lk_bdy )   CALL bdy_tra( kt )  ! BDY open boundaries
#endif
 
      ! set time step size (Euler/Leapfrog)
      IF( neuler == 0 .AND. kt == nit000 ) THEN   ;   r2dtra(:) =     rdttra(:)      ! at nit000             (Euler)
      ELSEIF( kt <= nit000 + 1 )           THEN   ;   r2dtra(:) = 2._wp* rdttra(:)      ! at nit000 or nit000+1 (Leapfrog)
      ENDIF

      ! trends computation initialisation
      IF( l_trdtra )   THEN                    ! store now fields before applying the Asselin filter
         CALL wrk_alloc( jpi, jpj, jpk, ztrdt, ztrds )
         ztrdt(:,:,:) = tsn(:,:,:,jp_tem) 
         ztrds(:,:,:) = tsn(:,:,:,jp_sal)
         IF( ln_traldf_iso ) THEN              ! diagnose the "pure" Kz diffusive trend 
            CALL trd_tra( kt, 'TRA', jp_tem, jptra_zdfp, ztrdt )
            CALL trd_tra( kt, 'TRA', jp_sal, jptra_zdfp, ztrds )
         ENDIF
      ENDIF

      IF( lk_diavecq ) THEN ! compute before and after thermal wind imbalance before Asselin filter
      IF ( l_diavecq_out)  THEN

        CALL wrk_alloc( jpi, jpj, jpk, zugb, zvgb, zuga, zvga)
        CALL wrk_alloc( jpi, jpj, jpk, zthub, zthvb, zthua, zthva, zdtthu, zdtthv )

        CALL geo_vel( tsb, zugb, zvgb, fsdept_n(:,:,:),kt) 
        CALL geo_vel( tsa, zuga, zvga, fsdept_n(:,:,:),kt) 
        CALL therm_wind_imb( ub_vecq, vb_vecq, zthub, zthvb, zugb, zvgb,kt)
        CALL therm_wind_imb( ua_vecq, va_vecq, zthua, zthva, zuga, zvga,kt)
        z1_2dt = 1._wp / (2. * rdt)
        zdtthu = ( zthua(:,:,:) - zthub(:,:,:) ) * z1_2dt
        zdtthv = ( zthva(:,:,:) - zthvb(:,:,:) ) * z1_2dt
        CALL iom_put( 'dtthu', zdtthu  )
        CALL iom_put( 'dtthv', zdtthv  )

        CALL wrk_dealloc( jpi, jpj, jpk, zugb, zvgb, zuga, zvga)
        CALL wrk_dealloc( jpi, jpj, jpk, zthub, zthvb, zthua, zthva, zdtthu, zdtthv )

      ENDIF
      ENDIF
 
      IF( neuler == 0 .AND. kt == nit000 ) THEN       ! Euler time-stepping at first time-step (only swap)
         DO jn = 1, jpts
            DO jk = 1, jpkm1
               tsn(:,:,jk,jn) = tsa(:,:,jk,jn)    
            END DO
         END DO
      ELSE                                            ! Leap-Frog + Asselin filter time stepping
         !
         IF( lk_vvl )  THEN   ;   CALL tra_nxt_vvl( kt, nit000, rdttra, 'TRA', tsb, tsn, tsa,   &
           &                                                              sbc_tsc, sbc_tsc_b, jpts )  ! variable volume level (vvl) 
         ELSE                 ;   CALL tra_nxt_fix( kt, nit000,         'TRA', tsb, tsn, tsa, jpts )  ! fixed    volume level 
         ENDIF
      ENDIF     
      !
     ! trends computation
      IF( l_trdtra ) THEN      ! trend of the Asselin filter (tb filtered - tb)/dt     
         DO jk = 1, jpkm1
            zfact = 1._wp / r2dtra(jk)             
            ztrdt(:,:,jk) = ( tsb(:,:,jk,jp_tem) - ztrdt(:,:,jk) ) * zfact
            ztrds(:,:,jk) = ( tsb(:,:,jk,jp_sal) - ztrds(:,:,jk) ) * zfact
         END DO
         CALL trd_tra( kt, 'TRA', jp_tem, jptra_atf, ztrdt )
         CALL trd_tra( kt, 'TRA', jp_sal, jptra_atf, ztrds )
         CALL wrk_dealloc( jpi, jpj, jpk, ztrdt, ztrds )
      END IF
      !
      !                        ! control print
      IF(ln_ctl)   CALL prt_ctl( tab3d_1=tsn(:,:,:,jp_tem), clinfo1=' nxt  - Tn: ', mask1=tmask,   &
         &                       tab3d_2=tsn(:,:,:,jp_sal), clinfo2=       ' Sn: ', mask2=tmask )
      !
      IF( nn_timing == 1 )   CALL timing_stop('tra_nxt')
      !
   END SUBROUTINE tra_nxt


   SUBROUTINE tra_nxt_fix( kt, kit000, cdtype, ptb, ptn, pta, kjpt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE tra_nxt_fix  ***
      !!
      !! ** Purpose :   fixed volume: apply the Asselin time filter and 
      !!                swap the tracer fields.
      !! 
      !! ** Method  : - Apply a Asselin time filter on now fields.
      !!              - save in (ta,sa) an average over the three time levels 
      !!             which will be used to compute rdn and thus the semi-implicit
      !!             hydrostatic pressure gradient (ln_dynhpg_imp = T)
      !!              - swap tracer fields to prepare the next time_step.
      !!                This can be summurized for tempearture as:
      !!             ztm = tn + rbcp * [ta -2 tn + tb ]       ln_dynhpg_imp = T
      !!             ztm = 0                                   otherwise
      !!                   with rbcp=1/4 * (1-atfp^4) / (1-atfp)
      !!             tb  = tn + atfp*[ tb - 2 tn + ta ]
      !!             tn  = ta  
      !!             ta  = ztm       (NB: reset to 0 after eos_bn2 call)
      !!
      !! ** Action  : - (tb,sb) and (tn,sn) ready for the next time step
      !!              - (ta,sa) time averaged (t,s)   (ln_dynhpg_imp = T)
      !!----------------------------------------------------------------------
      INTEGER         , INTENT(in   )                               ::   kt       ! ocean time-step index
      INTEGER         , INTENT(in   )                               ::   kit000   ! first time step index
      CHARACTER(len=3), INTENT(in   )                               ::   cdtype   ! =TRA or TRC (tracer indicator)
      INTEGER         , INTENT(in   )                               ::   kjpt     ! number of tracers
      REAL(wp)        , INTENT(inout), DIMENSION(jpi,jpj,jpk,kjpt)  ::   ptb      ! before tracer fields
      REAL(wp)        , INTENT(inout), DIMENSION(jpi,jpj,jpk,kjpt)  ::   ptn      ! now tracer fields
      REAL(wp)        , INTENT(inout), DIMENSION(jpi,jpj,jpk,kjpt)  ::   pta      ! tracer trend
      !
      INTEGER  ::   ji, jj, jk, jn   ! dummy loop indices
      LOGICAL  ::   ll_tra_hpg       ! local logical
      REAL(wp) ::   ztn, ztd         ! local scalars
      !!----------------------------------------------------------------------

      IF( kt == kit000 )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_nxt_fix : time stepping', cdtype
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
      ENDIF
      !
      IF( cdtype == 'TRA' )  THEN   ;   ll_tra_hpg = ln_dynhpg_imp    ! active  tracers case  and  semi-implicit hpg    
      ELSE                          ;   ll_tra_hpg = .FALSE.          ! passive tracers case or NO semi-implicit hpg
      ENDIF
      !
      DO jn = 1, kjpt
         !
         DO jk = 1, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  ztn = ptn(ji,jj,jk,jn)                                    
                  ztd = pta(ji,jj,jk,jn) - 2. * ztn + ptb(ji,jj,jk,jn)      !  time laplacian on tracers
                  !
                  ptb(ji,jj,jk,jn) = ztn + atfp * ztd                       ! ptb <-- filtered ptn 
                  ptn(ji,jj,jk,jn) = pta(ji,jj,jk,jn)                       ! ptn <-- pta
                  !
                  IF( ll_tra_hpg )   pta(ji,jj,jk,jn) = ztn + rbcp * ztd    ! pta <-- Brown & Campana average
               END DO
           END DO
         END DO
         !
      END DO
      !
   END SUBROUTINE tra_nxt_fix


   SUBROUTINE tra_nxt_vvl( kt, kit000, p2dt, cdtype, ptb, ptn, pta, psbc_tc, psbc_tc_b, kjpt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE tra_nxt_vvl  ***
      !!
      !! ** Purpose :   Time varying volume: apply the Asselin time filter  
      !!                and swap the tracer fields.
      !! 
      !! ** Method  : - Apply a thickness weighted Asselin time filter on now fields.
      !!              - save in (ta,sa) a thickness weighted average over the three 
      !!             time levels which will be used to compute rdn and thus the semi-
      !!             implicit hydrostatic pressure gradient (ln_dynhpg_imp = T)
      !!              - swap tracer fields to prepare the next time_step.
      !!                This can be summurized for tempearture as:
      !!             ztm = ( e3t_n*tn + rbcp*[ e3t_b*tb - 2 e3t_n*tn + e3t_a*ta ] )   ln_dynhpg_imp = T
      !!                  /( e3t_n    + rbcp*[ e3t_b    - 2 e3t_n    + e3t_a    ] )   
      !!             ztm = 0                                                       otherwise
      !!             tb  = ( e3t_n*tn + atfp*[ e3t_b*tb - 2 e3t_n*tn + e3t_a*ta ] )
      !!                  /( e3t_n    + atfp*[ e3t_b    - 2 e3t_n    + e3t_a    ] )
      !!             tn  = ta 
      !!             ta  = zt        (NB: reset to 0 after eos_bn2 call)
      !!
      !! ** Action  : - (tb,sb) and (tn,sn) ready for the next time step
      !!              - (ta,sa) time averaged (t,s)   (ln_dynhpg_imp = T)
      !!----------------------------------------------------------------------
      INTEGER         , INTENT(in   )                               ::  kt       ! ocean time-step index
      INTEGER         , INTENT(in   )                               ::  kit000   ! first time step index
      REAL(wp)        , INTENT(in   ), DIMENSION(jpk)               ::  p2dt     ! time-step
      CHARACTER(len=3), INTENT(in   )                               ::  cdtype   ! =TRA or TRC (tracer indicator)
      INTEGER         , INTENT(in   )                               ::  kjpt     ! number of tracers
      REAL(wp)        , INTENT(inout), DIMENSION(jpi,jpj,jpk,kjpt)  ::  ptb      ! before tracer fields
      REAL(wp)        , INTENT(inout), DIMENSION(jpi,jpj,jpk,kjpt)  ::  ptn      ! now tracer fields
      REAL(wp)        , INTENT(inout), DIMENSION(jpi,jpj,jpk,kjpt)  ::  pta      ! tracer trend
      REAL(wp)        , INTENT(in   ), DIMENSION(jpi,jpj,kjpt)      ::  psbc_tc   ! surface tracer content
      REAL(wp)        , INTENT(in   ), DIMENSION(jpi,jpj,kjpt)      ::  psbc_tc_b ! before surface tracer content

      !!     
      LOGICAL  ::   ll_tra_hpg, ll_traqsr, ll_rnf, ll_isf   ! local logical
      INTEGER  ::   ji, jj, jk, jn              ! dummy loop indices
      REAL(wp) ::   zfact1, ztc_a , ztc_n , ztc_b , ztc_f , ztc_d    ! local scalar
      REAL(wp) ::   zfact2, ze3t_b, ze3t_n, ze3t_a, ze3t_f, ze3t_d   !   -      -
      !!----------------------------------------------------------------------
      !
      IF( kt == kit000 )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_nxt_vvl : time stepping', cdtype
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
      ENDIF
      !
      IF( cdtype == 'TRA' )  THEN   
         ll_tra_hpg = ln_dynhpg_imp    ! active  tracers case  and  semi-implicit hpg
         ll_traqsr  = ln_traqsr        ! active  tracers case  and  solar penetration
         ll_rnf     = ln_rnf           ! active  tracers case  and  river runoffs
         IF (nn_isf .GE. 1) THEN 
            ll_isf = .TRUE.            ! active  tracers case  and  ice shelf melting/freezing
         ELSE
            ll_isf = .FALSE.
         END IF
      ELSE                          
         ll_tra_hpg = .FALSE.          ! passive tracers case or NO semi-implicit hpg
         ll_traqsr  = .FALSE.          ! active  tracers case and NO solar penetration
         ll_rnf     = .FALSE.          ! passive tracers or NO river runoffs
         ll_isf     = .FALSE.          ! passive tracers or NO ice shelf melting/freezing
      ENDIF
      !
      DO jn = 1, kjpt      
         DO jk = 1, jpkm1
            zfact1 = atfp * p2dt(jk)
            zfact2 = zfact1 / rau0
            DO jj = 1, jpj
               DO ji = 1, jpi
                  ze3t_b = fse3t_b(ji,jj,jk)
                  ze3t_n = fse3t_n(ji,jj,jk)
                  ze3t_a = fse3t_a(ji,jj,jk)
                  !                                         ! tracer content at Before, now and after
                  ztc_b  = ptb(ji,jj,jk,jn) * ze3t_b
                  ztc_n  = ptn(ji,jj,jk,jn) * ze3t_n
                  ztc_a  = pta(ji,jj,jk,jn) * ze3t_a
                  !
                  ze3t_d = ze3t_a - 2. * ze3t_n + ze3t_b
                  ztc_d  = ztc_a  - 2. * ztc_n  + ztc_b
                  !
                  ze3t_f = ze3t_n + atfp * ze3t_d
                  ztc_f  = ztc_n  + atfp * ztc_d
                  !
                  IF( jk == mikt(ji,jj) ) THEN           ! first level 
                     ze3t_f = ze3t_f - zfact2 * ( (emp_b(ji,jj)    - emp(ji,jj)   )  &
                            &                   - (rnf_b(ji,jj)    - rnf(ji,jj)   )  &
                            &                   + (fwfisf_b(ji,jj) - fwfisf(ji,jj))  )
                     ztc_f  = ztc_f  - zfact1 * ( psbc_tc(ji,jj,jn) - psbc_tc_b(ji,jj,jn) )
                  ENDIF

                  ! solar penetration (temperature only)
                  IF( ll_traqsr .AND. jn == jp_tem .AND. jk <= nksr )                            & 
                     &     ztc_f  = ztc_f  - zfact1 * ( qsr_hc(ji,jj,jk) - qsr_hc_b(ji,jj,jk) ) 

                  ! river runoff
                  IF( ll_rnf .AND. jk <= nk_rnf(ji,jj) )                                          &
                     &     ztc_f  = ztc_f  - zfact1 * ( rnf_tsc(ji,jj,jn) - rnf_tsc_b(ji,jj,jn) ) & 
                     &                              * fse3t_n(ji,jj,jk) / h_rnf(ji,jj)

                  ! ice shelf
                  IF( ll_isf ) THEN
                     ! level fully include in the Losch_2008 ice shelf boundary layer
                     IF ( jk >= misfkt(ji,jj) .AND. jk < misfkb(ji,jj) )                          &
                        ztc_f  = ztc_f  - zfact1 * ( risf_tsc(ji,jj,jn) - risf_tsc_b(ji,jj,jn) )  &
                               &                 * fse3t_n(ji,jj,jk) * r1_hisf_tbl (ji,jj)
                     ! level partially include in Losch_2008 ice shelf boundary layer 
                     IF ( jk == misfkb(ji,jj) )                                                   &
                        ztc_f  = ztc_f  - zfact1 * ( risf_tsc(ji,jj,jn) - risf_tsc_b(ji,jj,jn) )  &
                               &                 * fse3t_n(ji,jj,jk) * r1_hisf_tbl (ji,jj) * ralpha(ji,jj)
                  END IF

                  ze3t_f = 1.e0 / ze3t_f
                  ptb(ji,jj,jk,jn) = ztc_f * ze3t_f       ! ptb <-- ptn filtered
                  ptn(ji,jj,jk,jn) = pta(ji,jj,jk,jn)     ! ptn <-- pta
                  !
                  IF( ll_tra_hpg ) THEN        ! semi-implicit hpg (T & S only)
                     ze3t_d           = 1.e0   / ( ze3t_n + rbcp * ze3t_d )
                     pta(ji,jj,jk,jn) = ze3t_d * ( ztc_n  + rbcp * ztc_d  )   ! ta <-- Brown & Campana average
                  ENDIF
               END DO
            END DO
         END DO
         ! 
      END DO
      !
   END SUBROUTINE tra_nxt_vvl

   !!======================================================================
END MODULE tranxt
