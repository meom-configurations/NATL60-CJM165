MODULE trcldf
   !!======================================================================
   !!                       ***  MODULE  trcldf  ***
   !! Ocean Passive tracers : lateral diffusive trends
   !!=====================================================================
   !! History :  9.0  ! 2005-11 (G. Madec)  Original code
   !!       NEMO 3.0  ! 2008-01  (C. Ethe, G. Madec)  merge TRC-TRA
   !!----------------------------------------------------------------------
#if defined key_top
   !!----------------------------------------------------------------------
   !!   'key_top'                                                TOP models
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   trc_ldf     : update the tracer trend with the lateral diffusion
   !!       ldf_ctl : initialization, namelist read, and parameters control
   !!----------------------------------------------------------------------
   USE oce_trc         ! ocean dynamics and active tracers
   USE trc             ! ocean passive tracers variables
   USE trcnam_trp      ! passive tracers transport namelist variables
   USE ldftra_oce      ! lateral diffusion coefficient on tracers
   USE ldfslp          ! ???
   USE traldf_bilapg   ! lateral mixing            (tra_ldf_bilapg routine)
   USE traldf_bilap    ! lateral mixing            (tra_ldf_bilap routine)
   USE traldf_iso      ! lateral mixing            (tra_ldf_iso routine)
   USE traldf_iso_grif ! lateral mixing          (tra_ldf_iso_grif routine)
   USE traldf_lap      ! lateral mixing            (tra_ldf_lap routine)
   USE trd_oce
   USE trdtra
   USE prtctl_trc      ! Print control
   USE diasub
   USE trcdiasub


   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_ldf    ! called by step.F90
   !                                                 !!: ** lateral mixing namelist (nam_trcldf) **
   REAL(wp) ::  rldf_rat    ! ratio between active and passive tracers diffusive coefficient
   INTEGER  ::  nldf = 0   ! type of lateral diffusion used defined from ln_trcldf_... namlist logicals)
   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcldf.F90 6312 2016-02-15 11:43:52Z cetlod $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_ldf( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_ldf  ***
      !!
      !! ** Purpose :   compute the lateral ocean tracer physics.
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index
      !!
      INTEGER            :: ji, jj, jk, jn
      REAL(wp)           :: zdep
      CHARACTER (len=22) :: charout
      REAL(wp), POINTER, DIMENSION(:,:,:,:) ::   ztrtrd
      REAL(wp), POINTER, DIMENSION(:,:,:,:) ::   ztfu, ztfv, ztfw
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )   CALL timing_start('trc_ldf')
      !
      IF( kt == nittrc000 )   CALL ldf_ctl          ! initialisation & control of options

      rldf = rldf_rat
      !
      r_fact_lap(:,:,:) = 1.
      DO jk= 1, jpk
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( fsdept(ji,jj,jk) > 200. .AND. gphit(ji,jj) < 5. .AND. gphit(ji,jj) > -5. ) THEN
                  zdep = MAX( fsdept(ji,jj,jk) - 1000., 0. ) / 1000.
                  r_fact_lap(ji,jj,jk) = MAX( 1., rn_fact_lap * EXP( -zdep ) )
               ENDIF
            END DO
         END DO
      END DO
      !
      IF( l_trdtrc )  THEN
         CALL wrk_alloc( jpi, jpj, jpk, jptra, ztrtrd )
         ztrtrd(:,:,:,:)  = tra(:,:,:,:)
      ENDIF
      !
      IF( lk_diasub ) THEN
         CALL wrk_alloc( jpi, jpj, jpk, jptra, ztfu, ztfv, ztfw )
         ztfu(:,:,:,:)  = 0.   ;   ztfv(:,:,:,:)  = 0.   ;   ztfw(:,:,:,:)  = 0.
      ENDIF

      SELECT CASE ( nldf )                       ! compute lateral mixing trend and add it to the general trend
      CASE ( 0 )                                              ! iso-level laplacian
         IF( lk_diasub ) THEN ; CALL tra_ldf_lap( kt, nittrc000, 'TRC', gtru, gtrv, gtrui, gtrvi, trb, tra, jptra, ztfu, ztfv, ztfw )
         ELSE                 ; CALL tra_ldf_lap( kt, nittrc000, 'TRC', gtru, gtrv, gtrui, gtrvi, trb, tra, jptra )  
         ENDIF
      CASE ( 1 )     
         IF( ln_traldf_grif ) THEN
            CALL tra_ldf_iso_grif( kt, nittrc000, 'TRC', gtru, gtrv, trb, tra, jptra, rn_ahtb_0 )
         ELSE
            IF( lk_diasub ) THEN 
               CALL tra_ldf_iso( kt, nittrc000, 'TRC', gtru, gtrv, gtrui, gtrvi, trb, tra, jptra, rn_ahtb_0, ztfu, ztfv, ztfw ) 
            ELSE                 
               CALL tra_ldf_iso( kt, nittrc000, 'TRC', gtru, gtrv, gtrui, gtrvi, trb, tra, jptra, rn_ahtb_0 )
            ENDIF
         ENDIF
      CASE ( 2 )                                               ! iso-level bilaplacian
         IF( lk_diasub ) THEN ; CALL tra_ldf_bilap( kt, nittrc000, 'TRC', gtru, gtrv, gtrui, gtrvi, trb, tra, jptra, ztfu, ztfv, ztfw )
         ELSE                 ; CALL tra_ldf_bilap( kt, nittrc000, 'TRC', gtru, gtrv, gtrui, gtrvi, trb, tra, jptra )  
         ENDIF
      CASE ( 3 )   ;   CALL tra_ldf_bilapg( kt, nittrc000, 'TRC',             trb, tra, jptra            )  ! s-coord. horizontal bilaplacian
         !
      CASE ( -1 )                                     ! esopa: test all possibility with control print
         CALL tra_ldf_lap   ( kt, nittrc000, 'TRC', gtru, gtrv, gtrui, gtrvi, trb, tra, jptra            )
         WRITE(charout, FMT="('ldf0 ')") ;  CALL prt_ctl_trc_info(charout)
                                            CALL prt_ctl_trc( tab4d=tra, mask=tmask, clinfo=ctrcnm, clinfo2='trd' )
         IF( ln_traldf_grif ) THEN
            CALL tra_ldf_iso_grif( kt, nittrc000, 'TRC', gtru, gtrv, trb, tra, jptra, rn_ahtb_0 )
         ELSE
            CALL tra_ldf_iso     ( kt, nittrc000, 'TRC', gtru, gtrv, gtrui, gtrvi, trb, tra, jptra, rn_ahtb_0 )
         ENDIF
         WRITE(charout, FMT="('ldf1 ')") ;  CALL prt_ctl_trc_info(charout)
                                            CALL prt_ctl_trc( tab4d=tra, mask=tmask, clinfo=ctrcnm, clinfo2='trd' )
         CALL tra_ldf_bilap ( kt, nittrc000, 'TRC', gtru, gtrv, gtrui, gtrvi, trb, tra, jptra            )
         WRITE(charout, FMT="('ldf2 ')") ;  CALL prt_ctl_trc_info(charout)
                                            CALL prt_ctl_trc( tab4d=tra, mask=tmask, clinfo=ctrcnm, clinfo2='trd' )
         CALL tra_ldf_bilapg( kt, nittrc000, 'TRC',             trb, tra, jptra            )
         WRITE(charout, FMT="('ldf3 ')") ;  CALL prt_ctl_trc_info(charout)
                                            CALL prt_ctl_trc( tab4d=tra, mask=tmask, clinfo=ctrcnm, clinfo2='trd' )
      END SELECT
      !
      IF( l_trdtrc )   THEN                      ! save the horizontal diffusive trends for further diagnostics
        DO jn = 1, jptra
           ztrtrd(:,:,:,jn) = tra(:,:,:,jn) - ztrtrd(:,:,:,jn)
           CALL trd_tra( kt, 'TRC', jn, jptra_ldf, ztrtrd(:,:,:,jn) )
        END DO
        CALL wrk_dealloc( jpi, jpj, jpk, jptra, ztrtrd )
      ENDIF
      !                                          ! print mean trends (used for debugging)
      IF( lk_diasub ) THEN
        CALL trc_dia_sub_ldf( kt, ztfu, ztfv, ztfw )
        CALL wrk_dealloc( jpi, jpj, jpk, jptra, ztfu, ztfv, ztfw )
      ENDIF
      !                                               ! print mean trends (used for debugging)
      IF( ln_ctl )   THEN
         WRITE(charout, FMT="('ldf ')") ;  CALL prt_ctl_trc_info(charout)
                                           CALL prt_ctl_trc( tab4d=tra, mask=tmask, clinfo=ctrcnm, clinfo2='trd' )
      ENDIF
      !
      IF( nn_timing == 1 )   CALL timing_stop('trc_ldf')
      !
   END SUBROUTINE trc_ldf


   SUBROUTINE ldf_ctl
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ldf_ctl  ***
      !!
      !! ** Purpose :   Choice of the operator for the lateral tracer diffusion
      !!
      !! ** Method  :   set nldf from the namtra_ldf logicals
      !!      nldf == -2   No lateral diffusion
      !!      nldf == -1   ESOPA test: ALL operators are used
      !!      nldf ==  0   laplacian operator
      !!      nldf ==  1   Rotated laplacian operator
      !!      nldf ==  2   bilaplacian operator
      !!      nldf ==  3   Rotated bilaplacian
      !!----------------------------------------------------------------------
      INTEGER ::   ioptio, ierr         ! temporary integers
      !!----------------------------------------------------------------------

      IF (ABS(rn_aht_0) < 2._wp*TINY(1.e0)) THEN
         IF (ABS(rn_ahtrc_0) < 2._wp*TINY(1.e0)) THEN
            rldf_rat = 1.0_wp
         ELSE
            CALL ctl_stop( 'STOP', 'ldf_ctl : cannot define rldf_rat, rn_aht_0==0, rn_ahtrc_0 /=0' )
         END IF
      ELSE
         rldf_rat = rn_ahtrc_0 / rn_aht_0
      END IF
      !  Define the lateral mixing oparator for tracers
      ! ===============================================

      !                               ! control the input
      ioptio = 0
      IF( ln_trcldf_lap   )   ioptio = ioptio + 1
      IF( ln_trcldf_bilap )   ioptio = ioptio + 1
      IF( ioptio >  1 )   CALL ctl_stop( '          use ONE or NONE of the 2 lap/bilap operator type on tracer' )
      IF( ioptio == 0 )   nldf = -2   ! No lateral diffusion
      ioptio = 0
      IF( ln_trcldf_level )   ioptio = ioptio + 1
      IF( ln_trcldf_hor   )   ioptio = ioptio + 1
      IF( ln_trcldf_iso   )   ioptio = ioptio + 1
      IF( ioptio /= 1 )   CALL ctl_stop( '          use only ONE direction (level/hor/iso)' )

      ! defined the type of lateral diffusion from ln_trcldf_... logicals
      ! CAUTION : nldf = 1 is used in trazdf_imp, change it carefully
      ierr = 0
      IF( ln_trcldf_lap ) THEN       ! laplacian operator
         IF ( ln_zco ) THEN                ! z-coordinate
            IF ( ln_trcldf_level )   nldf = 0      ! iso-level  (no rotation)
            IF ( ln_trcldf_hor   )   nldf = 0      ! horizontal (no rotation)
            IF ( ln_trcldf_iso   )   nldf = 1      ! isoneutral (   rotation)
         ENDIF
         IF ( ln_zps ) THEN             ! z-coordinate
            IF ( ln_trcldf_level )   ierr = 1      ! iso-level not allowed
            IF ( ln_trcldf_hor   )   nldf = 0      ! horizontal (no rotation)
            IF ( ln_trcldf_iso   )   nldf = 1      ! isoneutral (   rotation)
         ENDIF
         IF ( ln_sco ) THEN             ! z-coordinate
            IF ( ln_trcldf_level )   nldf = 0      ! iso-level  (no rotation)
            IF ( ln_trcldf_hor   )   nldf = 1      ! horizontal (   rotation)
            IF ( ln_trcldf_iso   )   nldf = 1      ! isoneutral (   rotation)
         ENDIF
      ENDIF

      IF( ln_trcldf_bilap ) THEN      ! bilaplacian operator
         IF ( ln_zco ) THEN                ! z-coordinate
            IF ( ln_trcldf_level )   nldf = 2      ! iso-level  (no rotation)
            IF ( ln_trcldf_hor   )   nldf = 2      ! horizontal (no rotation)
            IF ( ln_trcldf_iso   )   ierr = 2      ! isoneutral (   rotation)
         ENDIF
         IF ( ln_zps ) THEN             ! z-coordinate
            IF ( ln_trcldf_level )   ierr = 1      ! iso-level not allowed
            IF ( ln_trcldf_hor   )   nldf = 2      ! horizontal (no rotation)
            IF ( ln_trcldf_iso   )   ierr = 2      ! isoneutral (   rotation)
         ENDIF
         IF ( ln_sco ) THEN             ! z-coordinate
            IF ( ln_trcldf_level )   nldf = 2      ! iso-level  (no rotation)
            IF ( ln_trcldf_hor   )   nldf = 3      ! horizontal (   rotation)
            IF ( ln_trcldf_iso   )   ierr = 2      ! isoneutral (   rotation)
         ENDIF
      ENDIF

      IF( ierr == 1 )   CALL ctl_stop( ' iso-level in z-coordinate - partial step, not allowed' )
      IF( ierr == 2 )   CALL ctl_stop( ' isoneutral bilaplacian operator does not exist' )
      IF( lk_traldf_eiv .AND. .NOT.ln_trcldf_iso )   &
           CALL ctl_stop( '          eddy induced velocity on tracers',   &
           &              ' the eddy induced velocity on tracers requires isopycnal laplacian diffusion' )
      IF( nldf == 1 .OR. nldf == 3 ) THEN      ! rotation
         IF( .NOT.lk_ldfslp )   CALL ctl_stop( '          the rotation of the diffusive tensor require key_ldfslp' )
#if defined key_offline
         l_traldf_rot = .TRUE.                 ! needed for trazdf_imp
#endif
      ENDIF

      IF( lk_esopa ) THEN
         IF(lwp) WRITE(numout,*) '          esopa control: use all lateral physics options'
         nldf = -1
      ENDIF

      IF(lwp) THEN
         WRITE(numout,*)
         IF( nldf == -2 )   WRITE(numout,*) '          NO lateral diffusion'
         IF( nldf == -1 )   WRITE(numout,*) '          ESOPA test All scheme used'
         IF( nldf ==  0 )   WRITE(numout,*) '          laplacian operator'
         IF( nldf ==  1 )   WRITE(numout,*) '          Rotated laplacian operator'
         IF( nldf ==  2 )   WRITE(numout,*) '          bilaplacian operator'
         IF( nldf ==  3 )   WRITE(numout,*) '          Rotated bilaplacian'
      ENDIF

      IF( ln_trcldf_bilap ) THEN
         IF(lwp) WRITE(numout,*) '          biharmonic tracer diffusion'
         IF( rn_ahtrc_0 > 0 .AND. .NOT. lk_esopa )   CALL ctl_stop( 'The horizontal diffusivity coef. rn_ahtrc_0 must be negative' )
      ELSE
         IF(lwp) WRITE(numout,*) '          harmonic tracer diffusion (default)'
         IF( rn_ahtrc_0 < 0 .AND. .NOT. lk_esopa )   CALL ctl_stop('The horizontal diffusivity coef. rn_ahtrc_0 must be positive' )
      ENDIF

      ! ratio between active and passive tracers diffusive coef.
      IF (ABS(rn_aht_0) < 2._wp*TINY(1.e0)) THEN
         IF (ABS(rn_ahtrc_0) < 2._wp*TINY(1.e0)) THEN
            rldf_rat = 1.0_wp
         ELSE
            CALL ctl_stop( 'STOP', 'ldf_ctl : cannot define rldf_rat, rn_aht_0==0, rn_ahtrc_0 /=0' )
         END IF
      ELSE
         rldf_rat = rn_ahtrc_0 / rn_aht_0
      END IF
      IF( rldf_rat < 0 ) THEN
         IF( .NOT.lk_offline ) THEN 
            CALL ctl_stop( 'Choose the same type of diffusive scheme both for active & passive tracers' )
         ELSE
            CALL ctl_stop( 'Change the sign of rn_aht_0 in namelist to -/+1' )
         ENDIF 
      ENDIF
      !
   END SUBROUTINE ldf_ctl
#else
   !!----------------------------------------------------------------------
   !!   Default option                                         Empty module
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_ldf( kt )
      INTEGER, INTENT(in) :: kt
      WRITE(*,*) 'trc_ldf: You should not have seen this print! error?', kt
   END SUBROUTINE trc_ldf
#endif
   !!======================================================================
END MODULE trcldf
