MODULE trcadv
   !!==============================================================================
   !!                       ***  MODULE  trcadv  ***
   !! Ocean passive tracers:  advection trend 
   !!==============================================================================
   !! History :  2.0  !  05-11  (G. Madec)  Original code
   !!            3.0  !  10-06  (C. Ethe)   Adapted to passive tracers
   !!----------------------------------------------------------------------
#if defined key_top
   !!----------------------------------------------------------------------
   !!   'key_top'                                                TOP models
   !!----------------------------------------------------------------------
   !!   trc_adv      : compute ocean tracer advection trend
   !!   trc_adv_ctl  : control the different options of advection scheme
   !!----------------------------------------------------------------------
   USE oce_trc         ! ocean dynamics and active tracers
   USE trc             ! ocean passive tracers variables
   USE trcnam_trp      ! passive tracers transport namelist variables
   USE traadv_cen2     ! 2nd order centered scheme (tra_adv_cen2   routine)
   USE traadv_tvd      ! TVD      scheme           (tra_adv_tvd    routine)
   USE traadv_muscl    ! MUSCL    scheme           (tra_adv_muscl  routine)
   USE traadv_muscl2   ! MUSCL2   scheme           (tra_adv_muscl2 routine)
   USE traadv_ubs      ! UBS      scheme           (tra_adv_ubs    routine)
   USE traadv_qck      ! QUICKEST scheme           (tra_adv_qck    routine)
   USE traadv_eiv      ! eddy induced velocity     (tra_adv_eiv    routine)
   USE traadv_mle      ! ML eddy induced velocity  (tra_adv_mle    routine)
   USE ldftra_oce      ! lateral diffusion coefficient on tracers
   USE prtctl_trc      ! Print control
   USE diasub
   USE trcdiasub


   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_adv          ! routine called by step module
   PUBLIC   trc_adv_alloc    ! routine called by nemogcm module

   INTEGER ::   nadv   ! choice of the type of advection scheme
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:) ::   r2dt  ! vertical profile time-step, = 2 rdttra
   !                                                    ! except at nitrrc000 (=rdttra) if neuler=0

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcadv.F90 5385 2015-06-09 13:50:42Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION trc_adv_alloc()
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trc_adv_alloc  ***
      !!----------------------------------------------------------------------

      ALLOCATE( r2dt(jpk), STAT=trc_adv_alloc )

      IF( trc_adv_alloc /= 0 ) CALL ctl_warn('trc_adv_alloc : failed to allocate array.')

   END FUNCTION trc_adv_alloc


   SUBROUTINE trc_adv( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trc_adv  ***
      !!
      !! ** Purpose :   compute the ocean tracer advection trend.
      !!
      !! ** Method  : - Update the tracer with the advection term following nadv
      !!----------------------------------------------------------------------
      !!
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      !
      INTEGER ::   jk 
      CHARACTER (len=22) ::   charout
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zun, zvn, zwn  ! effective velocity
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('trc_adv')
      !
      CALL wrk_alloc( jpi, jpj, jpk, zun, zvn, zwn )
      !

      IF( kt == nittrc000 )   CALL trc_adv_ctl          ! initialisation & control of options

      IF( ( neuler == 0 .AND. kt == nittrc000 ) .OR. ln_top_euler ) THEN     ! at nittrc000
         r2dt(:) =  rdttrc(:)           ! = rdttrc (use or restarting with Euler time stepping)
      ELSEIF( kt <= nittrc000 + nn_dttrc ) THEN          ! at nittrc000 or nittrc000+1
         r2dt(:) = 2. * rdttrc(:)       ! = 2 rdttrc (leapfrog)
      ENDIF
      !                                                   ! effective transport
      DO jk = 1, jpkm1
         !                                                ! eulerian transport only
         zun(:,:,jk) = e2u  (:,:) * fse3u(:,:,jk) * un(:,:,jk)
         zvn(:,:,jk) = e1v  (:,:) * fse3v(:,:,jk) * vn(:,:,jk)
         zwn(:,:,jk) = e1e2t(:,:)                 * wn(:,:,jk)
         !
      END DO
      !
      IF( ln_vvl_ztilde .OR. ln_vvl_layer ) THEN
         zun(:,:,:) = zun(:,:,:) + un_td(:,:,:)
         zvn(:,:,:) = zvn(:,:,:) + vn_td(:,:,:)
      ENDIF
      !
      zun(:,:,jpk) = 0._wp                                                     ! no transport trough the bottom
      zvn(:,:,jpk) = 0._wp                                                     ! no transport trough the bottom
      zwn(:,:,jpk) = 0._wp                                                     ! no transport trough the bottom


      IF( lk_diasub )  CALL trc_dia_sub_adv( kt, zun, zvn, zwn )  ! no eiv 

      IF( lk_traldf_eiv .AND. .NOT. ln_traldf_grif )  THEN  ! add the eiv transport (if necessary)
                          CALL tra_adv_eiv( kt, nittrc000, zun, zvn, zwn, 'TRC' )
         IF( lk_diasub )  CALL trc_dia_sub_adv_eiv( kt, zun, zvn, zwn )  
      ENDIF

      !
      IF( ln_mle    )   CALL tra_adv_mle( kt, nittrc000, zun, zvn, zwn, 'TRC' )    ! add the mle transport (if necessary)
      !
      SELECT CASE ( nadv )                            !==  compute advection trend and add it to general trend  ==!
      CASE ( 1 )   ;    CALL tra_adv_cen2  ( kt, nittrc000, 'TRC',       zun, zvn, zwn, trb, trn, tra, jptra )   !  2nd order centered
      CASE ( 2 )   ;    CALL tra_adv_tvd   ( kt, nittrc000, 'TRC', r2dt, zun, zvn, zwn, trb, trn, tra, jptra )   !  TVD 
      CASE ( 3 )   ;    CALL tra_adv_muscl ( kt, nittrc000, 'TRC', r2dt, zun, zvn, zwn, trb,      tra, jptra, ln_trcadv_msc_ups )   !  MUSCL 
      CASE ( 4 )   ;    CALL tra_adv_muscl2( kt, nittrc000, 'TRC', r2dt, zun, zvn, zwn, trb, trn, tra, jptra )   !  MUSCL2 
      CASE ( 5 )   ;    CALL tra_adv_ubs   ( kt, nittrc000, 'TRC', r2dt, zun, zvn, zwn, trb, trn, tra, jptra )   !  UBS 
      CASE ( 6 )   ;    CALL tra_adv_qck   ( kt, nittrc000, 'TRC', r2dt, zun, zvn, zwn, trb, trn, tra, jptra )   !  QUICKEST 
      !
      CASE (-1 )                                      !==  esopa: test all possibility with control print  ==!
         CALL tra_adv_cen2  ( kt, nittrc000, 'TRC',       zun, zvn, zwn, trb, trn, tra, jptra )          
         WRITE(charout, FMT="('adv1')")  ; CALL prt_ctl_trc_info(charout)
                                           CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm,clinfo2='trd')
         CALL tra_adv_tvd   ( kt, nittrc000, 'TRC', r2dt, zun, zvn, zwn, trb, trn, tra, jptra )          
         WRITE(charout, FMT="('adv2')")  ; CALL prt_ctl_trc_info(charout)
                                           CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm,clinfo2='trd')
         CALL tra_adv_muscl ( kt, nittrc000, 'TRC', r2dt, zun, zvn, zwn, trb,      tra, jptra, ln_trcadv_msc_ups  )          
         WRITE(charout, FMT="('adv3')")  ; CALL prt_ctl_trc_info(charout)
                                           CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm,clinfo2='trd')
         CALL tra_adv_muscl2( kt, nittrc000, 'TRC', r2dt, zun, zvn, zwn, trb, trn, tra, jptra )          
         WRITE(charout, FMT="('adv4')")  ; CALL prt_ctl_trc_info(charout)
                                           CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm,clinfo2='trd')
         CALL tra_adv_ubs   ( kt, nittrc000, 'TRC', r2dt, zun, zvn, zwn, trb, trn, tra, jptra )          
         WRITE(charout, FMT="('adv5')")  ; CALL prt_ctl_trc_info(charout)
                                           CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm,clinfo2='trd')
         CALL tra_adv_qck   ( kt, nittrc000, 'TRC', r2dt, zun, zvn, zwn, trb, trn, tra, jptra )          
         WRITE(charout, FMT="('adv6')")  ; CALL prt_ctl_trc_info(charout)
                                           CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm,clinfo2='trd')
         !
      END SELECT

      !                                              ! print mean trends (used for debugging)
      IF( ln_ctl )   THEN
         WRITE(charout, FMT="('adv ')")  ;  CALL prt_ctl_trc_info(charout)
                                            CALL prt_ctl_trc( tab4d=tra, mask=tmask, clinfo=ctrcnm, clinfo2='trd' )
      END IF
      !
      CALL wrk_dealloc( jpi, jpj, jpk, zun, zvn, zwn )
      !
      IF( nn_timing == 1 )  CALL timing_stop('trc_adv')
      !
   END SUBROUTINE trc_adv


   SUBROUTINE trc_adv_ctl
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE trc_adv_ctl  ***
      !!                
      !! ** Purpose : Control the consistency between namelist options for 
      !!              passive tracer advection schemes and set nadv
      !!----------------------------------------------------------------------
      INTEGER ::   ioptio
      !!----------------------------------------------------------------------
      IF(lwp) THEN                    ! Namelist print
         WRITE(numout,*)
         WRITE(numout,*) 'trc_adv_init : choice/control of the trccer advection scheme'
         WRITE(numout,*) '~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namtrc_adv : chose a advection scheme for trccers'
         WRITE(numout,*) '      2nd order advection scheme     ln_trcadv_cen2 = ', ln_trcadv_cen2
         WRITE(numout,*) '      TVD advection scheme           ln_trcadv_tvd = ', ln_trcadv_tvd
         WRITE(numout,*) '      MUSCL  advection scheme        ln_trcadv_muscl = ', ln_trcadv_muscl
         WRITE(numout,*) '      MUSCL2 advection scheme        ln_trcadv_muscl2 = ', ln_trcadv_muscl2
         WRITE(numout,*) '      UBS    advection scheme        ln_trcadv_ubs = ', ln_trcadv_ubs
         WRITE(numout,*) '      QUICKEST advection scheme      ln_trcadv_qck = ', ln_trcadv_qck
      ENDIF


      ioptio = 0                      ! Parameter control
      IF( ln_trcadv_cen2   )   ioptio = ioptio + 1
      IF( ln_trcadv_tvd    )   ioptio = ioptio + 1
      IF( ln_trcadv_muscl  )   ioptio = ioptio + 1
      IF( ln_trcadv_muscl2 )   ioptio = ioptio + 1
      IF( ln_trcadv_ubs    )   ioptio = ioptio + 1
      IF( ln_trcadv_qck    )   ioptio = ioptio + 1
      IF( lk_esopa         )   ioptio =          1

      IF( ioptio /= 1 )   CALL ctl_stop( 'Choose ONE advection scheme in namelist namtrc_adv' )

      !                              ! Set nadv
      IF( ln_trcadv_cen2   )   nadv =  1
      IF( ln_trcadv_tvd    )   nadv =  2
      IF( ln_trcadv_muscl  )   nadv =  3
      IF( ln_trcadv_muscl2 )   nadv =  4
      IF( ln_trcadv_ubs    )   nadv =  5
      IF( ln_trcadv_qck    )   nadv =  6
      IF( lk_esopa         )   nadv = -1

      IF(lwp) THEN                   ! Print the choice
         WRITE(numout,*)
         IF( nadv ==  1 )   WRITE(numout,*) '         2nd order scheme is used'
         IF( nadv ==  2 )   WRITE(numout,*) '         TVD       scheme is used'
         IF( nadv ==  3 )   WRITE(numout,*) '         MUSCL     scheme is used'
         IF( nadv ==  4 )   WRITE(numout,*) '         MUSCL2    scheme is used'
         IF( nadv ==  5 )   WRITE(numout,*) '         UBS       scheme is used'
         IF( nadv ==  6 )   WRITE(numout,*) '         QUICKEST  scheme is used'
         IF( nadv == -1 )   WRITE(numout,*) '         esopa test: use all advection scheme'
      ENDIF
      !
   END SUBROUTINE trc_adv_ctl
   
#else
   !!----------------------------------------------------------------------
   !!   Default option                                         Empty module
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_adv( kt )
      INTEGER, INTENT(in) :: kt
      WRITE(*,*) 'trc_adv: You should not have seen this print! error?', kt
   END SUBROUTINE trc_adv
#endif

  !!======================================================================
END MODULE trcadv
