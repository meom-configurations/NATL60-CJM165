MODULE diavecq
  !!======================================================================
  !!                       ***  MODULE  diavecq  ***
  !! vector Q components diagnostics
  !!======================================================================
  !! History :  3.2  !  2009-11  (S. Masson)  Original code
  !!            3.3  !  2010-10  (C. Ethe, G. Madec) reorganisation of initialisation phase + merge TRC-TRA
  !!            3.3  !  2016-05  (A. Albert) new diagnostics for vector Q
  !!----------------------------------------------------------------------
#if defined key_diavecq   
  !!----------------------------------------------------------------------
  !!   'key_diavecq'  :                       activate vector Q diagnotics
  !!----------------------------------------------------------------------
  !!   dia_vecq       : vector Q diagnostics
  !!   dia_vecq_init  : initialisation of vector Q diagnostics
  !!----------------------------------------------------------------------
  USE oce            ! ocean dynamics and active tracers 
  USE dom_oce        ! ocean space and time domain
  USE eosbn2         ! equation of state                (eos_bn2 routine)
  USE lib_mpp        ! distribued memory computing library
  USE iom            ! I/O manager library
  USE timing         ! preformance summary
  USE wrk_nemo       ! working arrays
  USE fldread        ! type FLD_N
  USE phycst         ! physical constant
  USE in_out_manager  ! I/O manager

  IMPLICIT NONE
  PRIVATE

  PUBLIC   dia_vecq        ! routine called in step.F90 module
  PUBLIC   geo_vel         ! routine called in tranxt.F90 module
  PUBLIC   therm_wind_imb  ! routine called in tranxt.F90 module

  LOGICAL, PUBLIC, PARAMETER :: lk_diavecq = .TRUE.   ! coupled flag
  LOGICAL, PUBLIC            :: l_diavecq_out = .FALSE.   ! coupled flag

  !! * Substitutions
#  include "domzgr_substitute.h90"
  !!----------------------------------------------------------------------
  !! NEMO/OPA 3.3 , NEMO Consortium (2010)
  !! $Id: diavecq.F90 5253 2015-05-05 12:46:46Z diovino $
  !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
  !!----------------------------------------------------------------------
CONTAINS

  SUBROUTINE dia_vecq( kt )
    !!----------------------------------------------------------------------
    !!                    ***  ROUTINE dia_vecq  ***
    !!
    !! ** Purpose :   compute and output vector Q components
    !!----------------------------------------------------------------------
    !
    INTEGER, INTENT( in ) ::   kt   ! ocean time-step index
    !
    INTEGER  ::   ji, jj, jk                      ! dummy loop arguments
    !
    REAL(wp) :: zarho, zmsv, zphv, zmsu, zphu, zalfg, zt
    REAL(wp), POINTER, DIMENSION(:,:,:)   :: zbuoy            ! 3D workspace
    REAL(wp), POINTER, DIMENSION(:,:,:)   :: zdxu             ! 3D workspace
    REAL(wp), POINTER, DIMENSION(:,:,:)   :: zdxv             ! 3D workspace
    REAL(wp), POINTER, DIMENSION(:,:,:)   :: zdyu             ! 3D workspace
    REAL(wp), POINTER, DIMENSION(:,:,:)   :: zdyv             ! 3D workspace
    REAL(wp), POINTER, DIMENSION(:,:,:)   :: zug              ! 3D workspace
    REAL(wp), POINTER, DIMENSION(:,:,:)   :: zvg              ! 3D workspace
    REAL(wp), POINTER, DIMENSION(:,:,:)   :: zthu             ! 3D workspace
    REAL(wp), POINTER, DIMENSION(:,:,:)   :: zthv             ! 3D workspace
    REAL(wp), POINTER, DIMENSION(:,:,:)   :: zdxug            ! 3D workspace
    REAL(wp), POINTER, DIMENSION(:,:,:)   :: zdxvg            ! 3D workspace
    REAL(wp), POINTER, DIMENSION(:,:,:)   :: zdyug            ! 3D workspace
    REAL(wp), POINTER, DIMENSION(:,:,:)   :: zdyvg            ! 3D workspace
    REAL(wp), POINTER, DIMENSION(:,:,:)   :: zrhd             ! 3D workspace
    REAL(wp), POINTER, DIMENSION(:,:,:)   :: zrhop            ! 3D workspace
    !!--------------------------------------------------------------------
!   IF( nn_timing == 1 )   CALL timing_start('dia_vecq')

    IF ( l_diavecq_out)  THEN

    CALL wrk_alloc( jpi , jpj, jpk         , zbuoy )
    CALL wrk_alloc( jpi , jpj, jpk         , zdxu )
    CALL wrk_alloc( jpi , jpj, jpk         , zdxv )
    CALL wrk_alloc( jpi , jpj, jpk         , zdyu )
    CALL wrk_alloc( jpi , jpj, jpk         , zdyv )
    CALL wrk_alloc( jpi , jpj, jpk         , zug )
    CALL wrk_alloc( jpi , jpj, jpk         , zvg )
    CALL wrk_alloc( jpi , jpj, jpk         , zthu )
    CALL wrk_alloc( jpi , jpj, jpk         , zthv )
    CALL wrk_alloc( jpi , jpj, jpk         , zdxug )
    CALL wrk_alloc( jpi , jpj, jpk         , zdxvg )
    CALL wrk_alloc( jpi , jpj, jpk         , zdyug )
    CALL wrk_alloc( jpi , jpj, jpk         , zdyvg )
    CALL wrk_alloc( jpi , jpj, jpk         , zrhd )
    CALL wrk_alloc( jpi , jpj, jpk         , zrhop )


       ! the buoyancy term
       CALL eos( tsn, zrhd, zrhop, fsdept_n(:,:,:) )
!      zbuoy(:,:,:) = 0._wp
       zbuoy(:,:,:) = -1 * grav / rau0 * zrhop(:,:,:)
!      CALL lbc_lnk( zbuoy, 'T', 1. )
       CALL iom_put( 'buoyancy', zbuoy )

       zt=kt
       ! horizontal gradient of current velocity
       zdxu(:,:,:) = 0._wp
       zdyu(:,:,:) = 0._wp
       zdxv(:,:,:) = 0._wp
       zdyv(:,:,:) = 0._wp

       DO jk = 1, jpk
          DO jj = 1, jpj
             DO ji = 2, jpim1
                zdxu(ji,jj,jk)= ( un(ji+1,jj,jk) - un(ji,jj,jk) ) / e1t(ji,jj)
                zdxv(ji,jj,jk)= ( vn(ji+1,jj,jk) - vn(ji,jj,jk) ) / e1f(ji,jj)
             END DO
          END DO

          DO jj = 2, jpjm1
             DO ji = 1, jpi
                zdyu(ji,jj,jk)= ( un(ji,jj+1,jk) - un(ji,jj,jk) ) / e2f(ji,jj)
                zdyv(ji,jj,jk)= ( vn(ji,jj+1,jk) - vn(ji,jj,jk) ) / e2t(ji,jj) 
             END DO
          END DO
       END DO

       CALL lbc_lnk( zdxu, 'T', 1. )
       CALL lbc_lnk( zdxv, 'F', 1. )
       CALL lbc_lnk( zdyu, 'F', 1. )
       CALL lbc_lnk( zdyv, 'T', 1. )

       CALL iom_put( 'dxu', zdxu )
       CALL iom_put( 'dxv', zdxv )
       CALL iom_put( 'dyu', zdyu )
       CALL iom_put( 'dyv', zdyv )


       ! horizontal  gradient of geostrophic currents

       !! geostrophic currents 
!      zug(:,:,:) = 0._wp
!      zvg(:,:,:) = 0._wp

!      CALL lbc_lnk( zug, 'U', -1. )
!      CALL lbc_lnk( zvg, 'V', -1. )

!   ENDIF !! l_diavecq_out 

    CALL geo_vel(tsn, zug, zvg, fsdept_n(:,:,:),kt)  ! always ? 

!   IF ( l_diavecq_out)  THEN
       !! horizontal gradient
       zdxug(:,:,:) = 0._wp
       zdyug(:,:,:) = 0._wp
       zdxvg(:,:,:) = 0._wp
       zdyvg(:,:,:) = 0._wp

       DO jk = 1, jpk
          DO jj = 1, jpj
             DO ji = 1, jpim1
                zdxug(ji,jj,jk)= ( zug(ji+1,jj,jk) - zug(ji,jj,jk) ) / e1t(ji,jj) 
                zdxvg(ji,jj,jk)= ( zvg(ji+1,jj,jk) - zvg(ji,jj,jk) ) / e1f(ji,jj) 
             END DO
          END DO

          DO jj = 1, jpjm1
             DO ji = 1, jpi
                zdyug(ji,jj,jk)= ( zug(ji,jj+1,jk) - zug(ji,jj,jk) ) / e2f(ji,jj) 
                zdyvg(ji,jj,jk)= ( zvg(ji,jj+1,jk) - zvg(ji,jj,jk) ) / e2t(ji,jj) 
             END DO
          END DO
       END DO

       CALL lbc_lnk( zdxug, 'T', 1. )
       CALL lbc_lnk( zdxvg, 'F', 1. )
       CALL lbc_lnk( zdyug, 'F', 1. )
       CALL lbc_lnk( zdyvg, 'T', 1. )

       CALL iom_put( 'dxug', zdxug )
       CALL iom_put( 'dxvg', zdxvg )
       CALL iom_put( 'dyug', zdyug )
       CALL iom_put( 'dyvg', zdyvg )

       ! thermal wind imbalance
!      zthu(:,:,:) = 0._wp
!      zthv(:,:,:) = 0._wp

!      CALL lbc_lnk( zthu, 'U', -1. )
!      CALL lbc_lnk( zthv, 'V', -1. )

!   ENDIF !! l_diavecq_out
    CALL therm_wind_imb( un, vn, zthu, zthv, zug, zvg, kt)

!   IF ( l_diavecq_out)  THEN

       CALL iom_put( 'thu', zthu )
       CALL iom_put( 'thv', zthv )

    !                     
    CALL wrk_dealloc( jpi , jpj, jpk         , zbuoy )
    CALL wrk_dealloc( jpi , jpj, jpk         , zdxu  )
    CALL wrk_dealloc( jpi , jpj, jpk         , zdxv  )
    CALL wrk_dealloc( jpi , jpj, jpk         , zdyu  )
    CALL wrk_dealloc( jpi , jpj, jpk         , zdyv  )
    CALL wrk_dealloc( jpi , jpj, jpk         , zug  )
    CALL wrk_dealloc( jpi , jpj, jpk         , zvg  )
    CALL wrk_dealloc( jpi , jpj, jpk         , zthu  )
    CALL wrk_dealloc( jpi , jpj, jpk         , zthv  )
    CALL wrk_dealloc( jpi , jpj, jpk         , zdxug  )
    CALL wrk_dealloc( jpi , jpj, jpk         , zdxvg  )
    CALL wrk_dealloc( jpi , jpj, jpk         , zdyug  )
    CALL wrk_dealloc( jpi , jpj, jpk         , zdyvg  )
    CALL wrk_dealloc( jpi , jpj, jpk         , zrhd  )
    CALL wrk_dealloc( jpi , jpj, jpk         , zrhop  )
    !
    ENDIF !! fin du if mod(kt-nit000+1,86400./rdt) = 0

!   IF( nn_timing == 1 )   CALL timing_stop('dia_vecq')
    !
  END SUBROUTINE dia_vecq

  SUBROUTINE geo_vel( pts, pug, pvg, pdep, kt)
    !!----------------------------------------------------------------------
    !!                   ***  ROUTINE geo_vel  ***
    !! ** Purpose :   Compute the i- and j- components of the geostrophic 
    !!       current in 3D
    !!
    !! ** Method  :  compute geostrophic component of the flow from
    !!               pressure gradient
    !!               results are respectively on U and V grid 
    !!----------------------------------------------------------------------
    REAL(wp), DIMENSION(jpi,jpj,jpk,jpts), INTENT(in   ) ::   pts   ! 1 : potential temperature  [Celcius]
    !                                                               ! 2 : salinity [psu]
    REAL(wp), DIMENSION(jpi,jpj,jpk     ), INTENT(  out) ::   pug   ! i- thermal wind imbalance
    REAL(wp), DIMENSION(jpi,jpj,jpk     ), INTENT(  out) ::   pvg   ! j- thermal wind imbalance
    REAL(wp), DIMENSION(jpi,jpj,jpk     ), INTENT(in   ) ::   pdep  ! depth [m]
    INTEGER                              , INTENT(in   ) ::   kt    ! ocean time-step index
    !
    INTEGER  ::   ji, jj, jk                ! dummy loop indices
    REAL(wp) ::   zt , zh , zs , ztm        ! local scalars
    REAL(wp) ::   zn , zn0, zn1, zn2, zn3   !   -      -
    REAL(wp) :: zarho, zmsv, zphv, zmsu, zphu, zalfg     
    REAL(wp), POINTER, DIMENSION(:,:,:)   :: zprn             ! 3D workspace
    REAL(wp), POINTER, DIMENSION(:,:,:)   :: zrhd             ! 3D workspace
    !!--------------------------------------------------------------------
!   IF( nn_timing == 1 )   CALL timing_start('geo_vel')

    CALL wrk_alloc( jpi , jpj, jpk         , zprn )
    CALL wrk_alloc( jpi , jpj, jpk         , zrhd )
    !!----------------------------------------------------------------------

!   IF ( l_diavecq_out)  THEN
       !! geostrophic currents 
       zprn(:,:,:) = 0._wp
       pug(:,:,:) = 0._wp
       pvg(:,:,:) = 0._wp

       zalfg = 0.5 * grav * rau0

       CALL eos( pts, zrhd, pdep )   ! ==> eos_insitu
       zprn(:,:,1) = zalfg * pdep(:,:,1) * ( 1 + zrhd(:,:,1) )       ! Surface value

       DO jk = 2, jpkm1                                              ! Vertical integration from the surface
          zprn(:,:,jk) = zprn(:,:,jk-1)   &
               &         + zalfg * pdep(:,:,jk) * ( 2. + zrhd(:,:,jk) + zrhd(:,:,jk-1) )
       END DO

       DO jk = 1, jpkm1
          DO jj = 2, jpjm1
             DO ji = 2, jpim1   ! vector opt.
                zmsv = 1. / MAX(  umask(ji-1,jj+1,jk) + umask(ji  ,jj+1,jk)   &
                                + umask(ji-1,jj  ,jk) + umask(ji  ,jj  ,jk) , 1. )
                zphv = ( zprn(ji  ,jj+1,jk) - zprn(ji-1,jj+1,jk) ) * umask(ji-1,jj+1,jk) / e1u(ji-1,jj+1)   &
                     + ( zprn(ji+1,jj+1,jk) - zprn(ji  ,jj+1,jk) ) * umask(ji  ,jj+1,jk) / e1u(ji  ,jj+1)   &
                     + ( zprn(ji  ,jj  ,jk) - zprn(ji-1,jj  ,jk) ) * umask(ji-1,jj  ,jk) / e1u(ji-1,jj  )   &
                     + ( zprn(ji+1,jj  ,jk) - zprn(ji  ,jj  ,jk) ) * umask(ji  ,jj  ,jk) / e1u(ji  ,jj  )
                zphv = 1. / rau0 * zphv * zmsv * vmask(ji,jj,jk)

                zmsu = 1. / MAX(  vmask(ji+1,jj  ,jk) + vmask(ji  ,jj  ,jk)   &
                                + vmask(ji+1,jj-1,jk) + vmask(ji  ,jj-1,jk) , 1. )
                zphu = ( zprn(ji+1,jj+1,jk) - zprn(ji+1,jj  ,jk) ) * vmask(ji+1,jj  ,jk) / e2v(ji+1,jj  )   &
                     + ( zprn(ji  ,jj+1,jk) - zprn(ji  ,jj  ,jk) ) * vmask(ji  ,jj  ,jk) / e2v(ji  ,jj  )   &
                     + ( zprn(ji+1,jj  ,jk) - zprn(ji+1,jj-1,jk) ) * vmask(ji+1,jj-1,jk) / e2v(ji+1,jj-1)   &
                     + ( zprn(ji  ,jj  ,jk) - zprn(ji  ,jj-1,jk) ) * vmask(ji  ,jj-1,jk) / e2v(ji  ,jj-1)
                zphu = 1. / rau0 * zphu * zmsu * umask(ji,jj,jk)

                pug(ji,jj,jk) = -2. * zphu / ( ff(ji,jj) + ff(ji  ,jj-1) )
                pvg(ji,jj,jk) =  2. * zphv / ( ff(ji,jj) + ff(ji-1,jj  ) )
             END DO
          END DO
       END DO

       CALL lbc_lnk( pug, 'U', -1. )
       CALL lbc_lnk( pvg, 'V', -1. )

!   ENDIF !! l_diavecq_out
    !                     
    CALL wrk_dealloc( jpi , jpj, jpk         , zprn  )
    CALL wrk_dealloc( jpi , jpj, jpk         , zrhd  )
    !                     
!   IF( nn_timing == 1 )   CALL timing_stop('geo_vel')
    !                     
  END SUBROUTINE geo_vel

  SUBROUTINE therm_wind_imb( pu, pv, pthu, pthv, pug, pvg, kt)
    !!----------------------------------------------------------------------
    !!                   ***  ROUTINE therm_wind_imb  ***
    !! ** Purpose :   Compute the i- and j- components of the thermal wind 
    !!       imbalance using u,v and geostrophic component
    !!
    !! ** Method  :   - compute horizontal gradient of the ageostrophic
    !!                  component times coriolis factor
    !!               results are respectively on UW and VW grid 
    !!----------------------------------------------------------------------
    REAL(wp), DIMENSION(jpi,jpj,jpk     ), INTENT(in   ) ::   pug   ! 
    REAL(wp), DIMENSION(jpi,jpj,jpk     ), INTENT(in   ) ::   pvg   ! 
    REAL(wp), DIMENSION(jpi,jpj,jpk     ), INTENT(in   ) ::   pu    ! i- thermal wind imbalance
    REAL(wp), DIMENSION(jpi,jpj,jpk     ), INTENT(in   ) ::   pv    ! j- thermal wind imbalance
    REAL(wp), DIMENSION(jpi,jpj,jpk     ), INTENT(  out) ::   pthu  ! i- thermal wind imbalance
    REAL(wp), DIMENSION(jpi,jpj,jpk     ), INTENT(  out) ::   pthv  ! j- thermal wind imbalance
    INTEGER                              , INTENT(in   ) ::   kt    ! ocean time-step index
    !
    INTEGER  ::   ji, jj, jk                ! dummy loop indices
    REAL(wp) :: zarho, zmsv, zphv, zmsu, zphu, zalfg     
    REAL(wp), POINTER, DIMENSION(:,:,:)   :: zua              ! 3D workspace
    REAL(wp), POINTER, DIMENSION(:,:,:)   :: zva              ! 3D workspace
    !!--------------------------------------------------------------------
!   IF( nn_timing == 1 )   CALL timing_start('therm_wind_imb')

    CALL wrk_alloc( jpi , jpj, jpk         , zua )
    CALL wrk_alloc( jpi , jpj, jpk         , zva )
    !!----------------------------------------------------------------------

!   IF ( l_diavecq_out)  THEN

       ! thermal wind imbalance
       zua(:,:,:) = 0._wp
       zva(:,:,:) = 0._wp
! JM tst
       pthu(:,:,:) = 0._wp
       pthv(:,:,:) = 0._wp

       DO jk = 1, jpkm1
          DO jj = 2, jpjm1
             DO ji = 2, jpim1  
                zua(ji,jj,jk) = pu(ji,jj,jk) - pug(ji,jj,jk)
                zva(ji,jj,jk) = pv(ji,jj,jk) - pvg(ji,jj,jk)
                pthu(ji,jj,jk) = ( ff(ji,jj) + ff(ji,jj-1) ) * ( zua(ji,jj,jk+1) - zua(ji,jj,jk) ) /  &
                        &        ( e3w_0(ji+1,jj,jk) + e3w_0(ji,jj,jk) )
                pthv(ji,jj,jk) = ( ff(ji,jj) + ff(ji-1,jj) ) * ( zva(ji,jj,jk+1) - zva(ji,jj,jk) ) /  &
                        &        ( e3w_0(ji,jj+1,jk) + e3w_0(ji,jj,jk) )
             END DO
          END DO
       END DO

       CALL lbc_lnk( pthu, 'UW', 1. ) 
       CALL lbc_lnk( pthv, 'VW', 1. ) 
       !                     
       CALL wrk_dealloc( jpi , jpj, jpk         , zua  )
       CALL wrk_dealloc( jpi , jpj, jpk         , zva  )
!   ENDIF !! l_diavecq_out
    !                     
!   IF( nn_timing == 1 )   CALL timing_stop('therm_wind_imb')
    !                     
  END SUBROUTINE therm_wind_imb

#else
  !!----------------------------------------------------------------------
  !!   Default option :                                         NO diavecq
  !!----------------------------------------------------------------------
  LOGICAL, PUBLIC, PARAMETER :: lk_diavecq = .FALSE.   ! coupled flag
  LOGICAL, PUBLIC            :: l_diavecq_out = .FALSE. ! coupled flag
CONTAINS
  SUBROUTINE dia_vecq_init    ! Dummy routine
  END SUBROUTINE dia_vecq_init
  SUBROUTINE dia_vecq( kt )   ! Empty routine
    INTEGER ::   kt
    WRITE(*,*) 'dia_vecq: You should not have seen this print! error?', kt
  END SUBROUTINE dia_vecq
#endif

  !!======================================================================
END MODULE diavecq
