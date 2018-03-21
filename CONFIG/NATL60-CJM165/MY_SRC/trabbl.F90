MODULE trabbl
   !!==============================================================================
   !!                       ***  MODULE  trabbl  ***
   !! Ocean physics :  advective and/or diffusive bottom boundary layer scheme
   !!==============================================================================
   !! History :  OPA  ! 1996-06  (L. Mortier)  Original code
   !!            8.0  ! 1997-11  (G. Madec)    Optimization
   !!   NEMO     1.0  ! 2002-08  (G. Madec)  free form + modules
   !!             -   ! 2004-01  (A. de Miranda, G. Madec, J.M. Molines ) add advective bbl
   !!            3.3  ! 2009-11  (G. Madec)  merge trabbl and trabbl_adv + style + optimization
   !!             -   ! 2010-04  (G. Madec)  Campin & Goosse advective bbl
   !!             -   ! 2010-06  (C. Ethe, G. Madec)  merge TRA-TRC
   !!             -   ! 2010-11  (G. Madec) add mbk. arrays associated to the deepest ocean level
   !!             -   ! 2013-04  (F. Roquet, G. Madec)  use of eosbn2 instead of local hard coded alpha and beta
   !!----------------------------------------------------------------------
#if   defined key_trabbl   ||   defined key_esopa
   !!----------------------------------------------------------------------
   !!   'key_trabbl'   or                             bottom boundary layer
   !!----------------------------------------------------------------------
   !!   tra_bbl_alloc : allocate trabbl arrays
   !!   tra_bbl       : update the tracer trends due to the bottom boundary layer (advective and/or diffusive)
   !!   tra_bbl_dif   : generic routine to compute bbl diffusive trend
   !!   tra_bbl_adv   : generic routine to compute bbl advective trend
   !!   bbl           : computation of bbl diffu. flux coef. & transport in bottom boundary layer
   !!   tra_bbl_init  : initialization, namelist read, parameters control
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and active tracers
   USE dom_oce        ! ocean space and time domain
   USE phycst         ! physical constant
   USE eosbn2         ! equation of state
   USE trd_oce     ! trends: ocean variables
   USE trdtra         ! trends: active tracers
   !
   USE iom            ! IOM library               
   USE in_out_manager ! I/O manager
   USE lbclnk         ! ocean lateral boundary conditions
   USE prtctl         ! Print control
   USE wrk_nemo       ! Memory Allocation
   USE timing         ! Timing
   USE lib_fortran    ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)  

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tra_bbl       !  routine called by step.F90
   PUBLIC   tra_bbl_init  !  routine called by opa.F90
   PUBLIC   tra_bbl_dif   !  routine called by trcbbl.F90
   PUBLIC   tra_bbl_adv   !  -          -          -              -
   PUBLIC   bbl           !  routine called by trcbbl.F90 and dtadyn.F90

   LOGICAL, PUBLIC, PARAMETER ::   lk_trabbl = .TRUE.    !: bottom boundary layer flag

   !                                !!* Namelist nambbl *
   INTEGER , PUBLIC ::   nn_bbl_ldf  !: =1   : diffusive bbl or not (=0)
   INTEGER , PUBLIC ::   nn_bbl_adv  !: =1/2 : advective bbl or not (=0)
   !                                            !  =1 : advective bbl using the bottom ocean velocity
   !                                            !  =2 :     -      -  using utr_bbl proportional to grad(rho)
   REAL(wp), PUBLIC ::   rn_ahtbbl   !: along slope bbl diffusive coefficient [m2/s]
   REAL(wp), PUBLIC ::   rn_gambbl   !: lateral coeff. for bottom boundary layer scheme [s]
!{ DRAKKAR
   LOGICAL , PUBLIC ::   ln_kriteria !: flag to use k-level instead of depth for bbl threshold
!}
   LOGICAL , PUBLIC ::   l_bbl                  !: flag to compute bbl diffu. flux coef and transport

   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:), PUBLIC ::   utr_bbl  , vtr_bbl   ! u- (v-) transport in the bottom boundary layer
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:), PUBLIC ::   ahu_bbl  , ahv_bbl   ! masked diffusive bbl coeff. at u & v-pts

   INTEGER , ALLOCATABLE, SAVE, DIMENSION(:,:), PUBLIC ::   mbku_d   , mbkv_d      ! vertical index of the "lower" bottom ocean U/V-level (PUBLIC for TAM)
   INTEGER , ALLOCATABLE, SAVE, DIMENSION(:,:), PUBLIC ::   mgrhu    , mgrhv       ! = +/-1, sign of grad(H) in u-(v-)direction (PUBLIC for TAM)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)         ::   ahu_bbl_0, ahv_bbl_0   ! diffusive bbl flux coefficients at u and v-points
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:), PUBLIC ::   e3u_bbl_0, e3v_bbl_0   ! thichness of the bbl (e3) at u and v-points (PUBLIC for TAM)

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: trabbl.F90 4990 2014-12-15 16:42:49Z timgraham $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION tra_bbl_alloc()
      !!----------------------------------------------------------------------
      !!                ***  FUNCTION tra_bbl_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( utr_bbl  (jpi,jpj) , ahu_bbl  (jpi,jpj) , mbku_d  (jpi,jpj) , mgrhu(jpi,jpj) ,     &
         &      vtr_bbl  (jpi,jpj) , ahv_bbl  (jpi,jpj) , mbkv_d  (jpi,jpj) , mgrhv(jpi,jpj) ,     &
         &      ahu_bbl_0(jpi,jpj) , ahv_bbl_0(jpi,jpj) ,                                          &
         &      e3u_bbl_0(jpi,jpj) , e3v_bbl_0(jpi,jpj) , STAT= tra_bbl_alloc                )
         !
      IF( lk_mpp            )   CALL mpp_sum ( tra_bbl_alloc )
      IF( tra_bbl_alloc > 0 )   CALL ctl_warn('tra_bbl_alloc: allocation of arrays failed.')
   END FUNCTION tra_bbl_alloc


   SUBROUTINE tra_bbl( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE bbl  ***
      !!
      !! ** Purpose :   Compute the before tracer (t & s) trend associated
      !!              with the bottom boundary layer and add it to the general
      !!              trend of tracer equations.
      !!
      !! ** Method  :   Depending on namtra_bbl namelist parameters the bbl
      !!              diffusive and/or advective contribution to the tracer trend
      !!              is added to the general tracer trend
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step
      !
      REAL(wp), POINTER, DIMENSION(:,:,:) ::  ztrdt, ztrds
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start( 'tra_bbl')
      !
      IF( l_trdtra )   THEN                         !* Save ta and sa trends
         CALL wrk_alloc( jpi, jpj, jpk, ztrdt, ztrds )
         ztrdt(:,:,:) = tsa(:,:,:,jp_tem)
         ztrds(:,:,:) = tsa(:,:,:,jp_sal)
      ENDIF

      IF( l_bbl )   CALL bbl( kt, nit000, 'TRA' )   !* bbl coef. and transport (only if not already done in trcbbl)

      IF( nn_bbl_ldf == 1 ) THEN                    !* Diffusive bbl
         !
         CALL tra_bbl_dif( tsb, tsa, jpts )
         IF( ln_ctl )  &
         CALL prt_ctl( tab3d_1=tsa(:,:,:,jp_tem), clinfo1=' bbl_ldf  - Ta: ', mask1=tmask, &
            &          tab3d_2=tsa(:,:,:,jp_sal), clinfo2=           ' Sa: ', mask2=tmask, clinfo3='tra' )
         ! lateral boundary conditions ; just need for outputs
         CALL lbc_lnk( ahu_bbl, 'U', 1. )     ;     CALL lbc_lnk( ahv_bbl, 'V', 1. )
         CALL iom_put( "ahu_bbl", ahu_bbl )   ! bbl diffusive flux i-coef
         CALL iom_put( "ahv_bbl", ahv_bbl )   ! bbl diffusive flux j-coef
         !
      END IF

      IF( nn_bbl_adv /= 0 ) THEN                    !* Advective bbl
         !
         CALL tra_bbl_adv( tsb, tsa, jpts )
         IF(ln_ctl)   &
         CALL prt_ctl( tab3d_1=tsa(:,:,:,jp_tem), clinfo1=' bbl_adv  - Ta: ', mask1=tmask,   &
            &          tab3d_2=tsa(:,:,:,jp_sal), clinfo2=           ' Sa: ', mask2=tmask, clinfo3='tra' )
         ! lateral boundary conditions ; just need for outputs
         CALL lbc_lnk( utr_bbl, 'U', 1. )     ;   CALL lbc_lnk( vtr_bbl, 'V', 1. )
         CALL iom_put( "uoce_bbl", utr_bbl )  ! bbl i-transport
         CALL iom_put( "voce_bbl", vtr_bbl )  ! bbl j-transport
         !
      END IF

      IF( l_trdtra )   THEN                      ! save the horizontal diffusive trends for further diagnostics
         ztrdt(:,:,:) = tsa(:,:,:,jp_tem) - ztrdt(:,:,:)
         ztrds(:,:,:) = tsa(:,:,:,jp_sal) - ztrds(:,:,:)
         CALL trd_tra( kt, 'TRA', jp_tem, jptra_bbl, ztrdt )
         CALL trd_tra( kt, 'TRA', jp_sal, jptra_bbl, ztrds )
         CALL wrk_dealloc( jpi, jpj, jpk, ztrdt, ztrds )
      ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop( 'tra_bbl')
      !
   END SUBROUTINE tra_bbl


   SUBROUTINE tra_bbl_dif( ptb, pta, kjpt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_bbl_dif  ***
      !!
      !! ** Purpose :   Computes the bottom boundary horizontal and vertical
      !!                advection terms.
      !!
      !! ** Method  : * diffusive bbl only (nn_bbl_ldf=1) :
      !!        When the product grad( rho) * grad(h) < 0 (where grad is an
      !!      along bottom slope gradient) an additional lateral 2nd order
      !!      diffusion along the bottom slope is added to the general
      !!      tracer trend, otherwise the additional trend is set to 0.
      !!      A typical value of ahbt is 2000 m2/s (equivalent to
      !!      a downslope velocity of 20 cm/s if the condition for slope
      !!      convection is satified)
      !!
      !! ** Action  :   pta   increased by the bbl diffusive trend
      !!
      !! References : Beckmann, A., and R. Doscher, 1997, J. Phys.Oceanogr., 581-591.
      !!              Campin, J.-M., and H. Goosse, 1999, Tellus, 412-430.
      !!----------------------------------------------------------------------
      INTEGER                              , INTENT(in   ) ::   kjpt   ! number of tracers
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(in   ) ::   ptb    ! before and now tracer fields
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(inout) ::   pta    ! tracer trend
      !
      INTEGER  ::   ji, jj, jn   ! dummy loop indices
      INTEGER  ::   ik           ! local integers
      REAL(wp) ::   zbtr         ! local scalars
      REAL(wp), POINTER, DIMENSION(:,:) :: zptb
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('tra_bbl_dif')
      !
      CALL wrk_alloc( jpi, jpj, zptb )
      !
      DO jn = 1, kjpt                                     ! tracer loop
         !                                                ! ===========
         DO jj = 1, jpj
            DO ji = 1, jpi
               ik = mbkt(ji,jj)                              ! bottom T-level index
               zptb(ji,jj) = ptb(ji,jj,ik,jn)       ! bottom before T and S
            END DO
         END DO
         !               
         DO jj = 2, jpjm1                                    ! Compute the trend
            DO ji = 2, jpim1
               ik = mbkt(ji,jj)                              ! bottom T-level index
               zbtr = r1_e12t(ji,jj)  / fse3t(ji,jj,ik)
               pta(ji,jj,ik,jn) = pta(ji,jj,ik,jn)                                                         &
                  &               + (   ahu_bbl(ji  ,jj  ) * ( zptb(ji+1,jj  ) - zptb(ji  ,jj  ) )   &
                  &                   - ahu_bbl(ji-1,jj  ) * ( zptb(ji  ,jj  ) - zptb(ji-1,jj  ) )   &
                  &                   + ahv_bbl(ji  ,jj  ) * ( zptb(ji  ,jj+1) - zptb(ji  ,jj  ) )   &
                  &                   - ahv_bbl(ji  ,jj-1) * ( zptb(ji  ,jj  ) - zptb(ji  ,jj-1) )   ) * zbtr
            END DO
         END DO
         !                                                  ! ===========
      END DO                                                ! end tracer
      !                                                     ! ===========
      CALL wrk_dealloc( jpi, jpj, zptb )
      !
      IF( nn_timing == 1 )  CALL timing_stop('tra_bbl_dif')
      !
   END SUBROUTINE tra_bbl_dif


   SUBROUTINE tra_bbl_adv( ptb, pta, kjpt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trc_bbl  ***
      !!
      !! ** Purpose :   Compute the before passive tracer trend associated
      !!     with the bottom boundary layer and add it to the general trend
      !!     of tracer equations.
      !! ** Method  :   advective bbl (nn_bbl_adv = 1 or 2) :
      !!      nn_bbl_adv = 1   use of the ocean near bottom velocity as bbl velocity
      !!      nn_bbl_adv = 2   follow Campin and Goosse (1999) implentation i.e.
      !!                       transport proportional to the along-slope density gradient
      !!
      !! References : Beckmann, A., and R. Doscher, 1997, J. Phys.Oceanogr., 581-591.
      !!              Campin, J.-M., and H. Goosse, 1999, Tellus, 412-430.
      !!----------------------------------------------------------------------
      INTEGER                              , INTENT(in   ) ::   kjpt   ! number of tracers
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(in   ) ::   ptb    ! before and now tracer fields
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(inout) ::   pta    ! tracer trend
      !
      INTEGER  ::   ji, jj, jk, jn           ! dummy loop indices
      INTEGER  ::   iis , iid , ijs , ijd    ! local integers
      INTEGER  ::   ikus, ikud, ikvs, ikvd   !   -       -
      REAL(wp) ::   zbtr, ztra               ! local scalars
      REAL(wp) ::   zu_bbl, zv_bbl           !   -      -
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start( 'tra_bbl_adv')
      !                                                          ! ===========
      DO jn = 1, kjpt                                            ! tracer loop
         !                                                       ! ===========
         DO jj = 1, jpjm1
            DO ji = 1, jpim1            ! CAUTION start from i=1 to update i=2 when cyclic east-west
               IF( utr_bbl(ji,jj) /= 0.e0 ) THEN            ! non-zero i-direction bbl advection
                  ! down-slope i/k-indices (deep)      &   up-slope i/k indices (shelf)
                  iid  = ji + MAX( 0, mgrhu(ji,jj) )   ;   iis  = ji + 1 - MAX( 0, mgrhu(ji,jj) )
                  ikud = mbku_d(ji,jj)                 ;   ikus = mbku(ji,jj)
                  zu_bbl = ABS( utr_bbl(ji,jj) )
                  !
                  !                                               ! up  -slope T-point (shelf bottom point)
                  zbtr = r1_e12t(iis,jj) / fse3t(iis,jj,ikus)
                  ztra = zu_bbl * ( ptb(iid,jj,ikus,jn) - ptb(iis,jj,ikus,jn) ) * zbtr
                  pta(iis,jj,ikus,jn) = pta(iis,jj,ikus,jn) + ztra
                  !
                  DO jk = ikus, ikud-1                            ! down-slope upper to down T-point (deep column)
                     zbtr = r1_e12t(iid,jj) / fse3t(iid,jj,jk)
                     ztra = zu_bbl * ( ptb(iid,jj,jk+1,jn) - ptb(iid,jj,jk,jn) ) * zbtr
                     pta(iid,jj,jk,jn) = pta(iid,jj,jk,jn) + ztra
                  END DO
                  !
                  zbtr = r1_e12t(iid,jj) / fse3t(iid,jj,ikud)
                  ztra = zu_bbl * ( ptb(iis,jj,ikus,jn) - ptb(iid,jj,ikud,jn) ) * zbtr
                  pta(iid,jj,ikud,jn) = pta(iid,jj,ikud,jn) + ztra
               ENDIF
               !
               IF( vtr_bbl(ji,jj) /= 0.e0 ) THEN            ! non-zero j-direction bbl advection
                  ! down-slope j/k-indices (deep)        &   up-slope j/k indices (shelf)
                  ijd  = jj + MAX( 0, mgrhv(ji,jj) )     ;   ijs  = jj + 1 - MAX( 0, mgrhv(ji,jj) )
                  ikvd = mbkv_d(ji,jj)                   ;   ikvs = mbkv(ji,jj)
                  zv_bbl = ABS( vtr_bbl(ji,jj) )
                  !
                  ! up  -slope T-point (shelf bottom point)
                  zbtr = r1_e12t(ji,ijs) / fse3t(ji,ijs,ikvs)
                  ztra = zv_bbl * ( ptb(ji,ijd,ikvs,jn) - ptb(ji,ijs,ikvs,jn) ) * zbtr
                  pta(ji,ijs,ikvs,jn) = pta(ji,ijs,ikvs,jn) + ztra
                  !
                  DO jk = ikvs, ikvd-1                            ! down-slope upper to down T-point (deep column)
                     zbtr = r1_e12t(ji,ijd) / fse3t(ji,ijd,jk)
                     ztra = zv_bbl * ( ptb(ji,ijd,jk+1,jn) - ptb(ji,ijd,jk,jn) ) * zbtr
                     pta(ji,ijd,jk,jn) = pta(ji,ijd,jk,jn)  + ztra
                  END DO
                  !                                               ! down-slope T-point (deep bottom point)
                  zbtr = r1_e12t(ji,ijd) / fse3t(ji,ijd,ikvd)
                  ztra = zv_bbl * ( ptb(ji,ijs,ikvs,jn) - ptb(ji,ijd,ikvd,jn) ) * zbtr
                  pta(ji,ijd,ikvd,jn) = pta(ji,ijd,ikvd,jn) + ztra
               ENDIF
            END DO
            !
         END DO
         !                                                  ! ===========
      END DO                                                ! end tracer
      !                                                     ! ===========
      !
      IF( nn_timing == 1 )  CALL timing_stop( 'tra_bbl_adv')
      !
   END SUBROUTINE tra_bbl_adv


   SUBROUTINE bbl( kt, kit000, cdtype )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE bbl  ***
      !!
      !! ** Purpose :   Computes the bottom boundary horizontal and vertical
      !!                advection terms.
      !!
      !! ** Method  : * diffusive bbl (nn_bbl_ldf=1) :
      !!        When the product grad( rho) * grad(h) < 0 (where grad is an
      !!      along bottom slope gradient) an additional lateral 2nd order
      !!      diffusion along the bottom slope is added to the general
      !!      tracer trend, otherwise the additional trend is set to 0.
      !!      A typical value of ahbt is 2000 m2/s (equivalent to
      !!      a downslope velocity of 20 cm/s if the condition for slope
      !!      convection is satified)
      !!              * advective bbl (nn_bbl_adv=1 or 2) :
      !!      nn_bbl_adv = 1   use of the ocean velocity as bbl velocity
      !!      nn_bbl_adv = 2   follow Campin and Goosse (1999) implentation
      !!        i.e. transport proportional to the along-slope density gradient
      !!
      !!      NB: the along slope density gradient is evaluated using the
      !!      local density (i.e. referenced at a common local depth).
      !!
      !! References : Beckmann, A., and R. Doscher, 1997, J. Phys.Oceanogr., 581-591.
      !!              Campin, J.-M., and H. Goosse, 1999, Tellus, 412-430.
      !!----------------------------------------------------------------------
      INTEGER         , INTENT(in   ) ::   kt       ! ocean time-step index
      INTEGER         , INTENT(in   ) ::   kit000   ! first time step index
      CHARACTER(len=3), INTENT(in   ) ::   cdtype   ! =TRA or TRC (tracer indicator)
      !!
      INTEGER  ::   ji, jj                    ! dummy loop indices
      INTEGER  ::   ik                        ! local integers
      INTEGER  ::   iis, iid, ikus, ikud      !   -       -
      INTEGER  ::   ijs, ijd, ikvs, ikvd      !   -       -
      REAL(wp) ::   za, zb, zgdrho            ! local scalars
      REAL(wp) ::   zsign, zsigna, zgbbl      !   -      -
      REAL(wp), DIMENSION(jpi,jpj,jpts)   :: zts, zab         ! 3D workspace
      REAL(wp), DIMENSION(jpi,jpj)        :: zub, zvb, zdep   ! 2D workspace
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start( 'bbl')
      !
      IF( kt == kit000 )  THEN
         IF(lwp)  WRITE(numout,*)
         IF(lwp)  WRITE(numout,*) 'trabbl:bbl : Compute bbl velocities and diffusive coefficients in ', cdtype
         IF(lwp)  WRITE(numout,*) '~~~~~~~~~~'
      ENDIF
      !                                        !* bottom variables (T, S, alpha, beta, depth, velocity)
      DO jj = 1, jpj
         DO ji = 1, jpi
            ik = mbkt(ji,jj)                             ! bottom T-level index
            zts (ji,jj,jp_tem) = tsb(ji,jj,ik,jp_tem)    ! bottom before T and S
            zts (ji,jj,jp_sal) = tsb(ji,jj,ik,jp_sal)
            !
            zdep(ji,jj) = fsdept(ji,jj,ik)               ! bottom T-level reference depth
            zub (ji,jj) = un(ji,jj,mbku(ji,jj))          ! bottom velocity
            zvb (ji,jj) = vn(ji,jj,mbkv(ji,jj))
         END DO
      END DO
      !
      CALL eos_rab( zts, zdep, zab )
      !
      !                                   !-------------------!
      IF( nn_bbl_ldf == 1 ) THEN          !   diffusive bbl   !
         !                                !-------------------!
         DO jj = 1, jpjm1                      ! (criteria for non zero flux: grad(rho).grad(h) < 0 )
            DO ji = 1, fs_jpim1   ! vector opt.
               !                                                   ! i-direction
               za = zab(ji+1,jj,jp_tem) + zab(ji,jj,jp_tem)              ! 2*(alpha,beta) at u-point
               zb = zab(ji+1,jj,jp_sal) + zab(ji,jj,jp_sal)
               !                                                         ! 2*masked bottom density gradient
               zgdrho = (  za * ( zts(ji+1,jj,jp_tem) - zts(ji,jj,jp_tem) )    &
                  &      - zb * ( zts(ji+1,jj,jp_sal) - zts(ji,jj,jp_sal) )  ) * umask(ji,jj,1)
               !
               zsign  = SIGN(  0.5, -zgdrho * REAL( mgrhu(ji,jj) )  )    ! sign of ( i-gradient * i-slope )
               ahu_bbl(ji,jj) = ( 0.5 - zsign ) * ahu_bbl_0(ji,jj)       ! masked diffusive flux coeff.
               !
               !                                                   ! j-direction
               za = zab(ji,jj+1,jp_tem) + zab(ji,jj,jp_tem)              ! 2*(alpha,beta) at v-point
               zb = zab(ji,jj+1,jp_sal) + zab(ji,jj,jp_sal)
               !                                                         ! 2*masked bottom density gradient
               zgdrho = (  za * ( zts(ji,jj+1,jp_tem) - zts(ji,jj,jp_tem) )    &
                  &      - zb * ( zts(ji,jj+1,jp_sal) - zts(ji,jj,jp_sal) )  ) * vmask(ji,jj,1)
               !
               zsign = SIGN(  0.5, -zgdrho * REAL( mgrhv(ji,jj) )  )     ! sign of ( j-gradient * j-slope )
               ahv_bbl(ji,jj) = ( 0.5 - zsign ) * ahv_bbl_0(ji,jj)
            END DO
         END DO
         !
      ENDIF

      !                                   !-------------------!
      IF( nn_bbl_adv /= 0 ) THEN          !   advective bbl   !
         !                                !-------------------!
         SELECT CASE ( nn_bbl_adv )             !* bbl transport type
         !
         CASE( 1 )                                   != use of upper velocity
            DO jj = 1, jpjm1                                 ! criteria: grad(rho).grad(h)<0  and grad(rho).grad(h)<0
               DO ji = 1, fs_jpim1   ! vector opt.
                  !                                                  ! i-direction
                  za = zab(ji+1,jj,jp_tem) + zab(ji,jj,jp_tem)               ! 2*(alpha,beta) at u-point
                  zb = zab(ji+1,jj,jp_sal) + zab(ji,jj,jp_sal)
                  !                                                          ! 2*masked bottom density gradient 
                  zgdrho = (  za * ( zts(ji+1,jj,jp_tem) - zts(ji,jj,jp_tem) )    &
                            - zb * ( zts(ji+1,jj,jp_sal) - zts(ji,jj,jp_sal) )  ) * umask(ji,jj,1)
                  !
                  zsign = SIGN(  0.5, - zgdrho   * REAL( mgrhu(ji,jj) )  )   ! sign of i-gradient * i-slope
                  zsigna= SIGN(  0.5, zub(ji,jj) * REAL( mgrhu(ji,jj) )  )   ! sign of u * i-slope
                  !
                  !                                                          ! bbl velocity
                  utr_bbl(ji,jj) = ( 0.5 + zsigna ) * ( 0.5 - zsign ) * e2u(ji,jj) * e3u_bbl_0(ji,jj) * zub(ji,jj)
                  !
                  !                                                  ! j-direction
                  za = zab(ji,jj+1,jp_tem) + zab(ji,jj,jp_tem)               ! 2*(alpha,beta) at v-point
                  zb = zab(ji,jj+1,jp_sal) + zab(ji,jj,jp_sal)
                  !                                                          ! 2*masked bottom density gradient
                  zgdrho = (  za * ( zts(ji,jj+1,jp_tem) - zts(ji,jj,jp_tem) )    &
                     &      - zb * ( zts(ji,jj+1,jp_sal) - zts(ji,jj,jp_sal) )  ) * vmask(ji,jj,1)
                  zsign = SIGN(  0.5, - zgdrho   * REAL( mgrhv(ji,jj) )  )   ! sign of j-gradient * j-slope
                  zsigna= SIGN(  0.5, zvb(ji,jj) * REAL( mgrhv(ji,jj) )  )   ! sign of u * i-slope
                  !
                  !                                                          ! bbl transport
                  vtr_bbl(ji,jj) = ( 0.5 + zsigna ) * ( 0.5 - zsign ) * e1v(ji,jj) * e3v_bbl_0(ji,jj) * zvb(ji,jj)
               END DO
            END DO
            !
         CASE( 2 )                                 != bbl velocity = F( delta rho )
            zgbbl = grav * rn_gambbl
            DO jj = 1, jpjm1                            ! criteria: rho_up > rho_down
               DO ji = 1, fs_jpim1   ! vector opt.
                  !                                                  ! i-direction
                  ! down-slope T-point i/k-index (deep)  &   up-slope T-point i/k-index (shelf)
                  iid  = ji + MAX( 0, mgrhu(ji,jj) )
                  iis  = ji + 1 - MAX( 0, mgrhu(ji,jj) )
                  !
                  ikud = mbku_d(ji,jj)
                  ikus = mbku(ji,jj)
                  !
                  za = zab(ji+1,jj,jp_tem) + zab(ji,jj,jp_tem)               ! 2*(alpha,beta) at u-point
                  zb = zab(ji+1,jj,jp_sal) + zab(ji,jj,jp_sal)
                  !                                                          !   masked bottom density gradient
                  zgdrho = 0.5 * (  za * ( zts(iid,jj,jp_tem) - zts(iis,jj,jp_tem) )    &
                     &            - zb * ( zts(iid,jj,jp_sal) - zts(iis,jj,jp_sal) )  ) * umask(ji,jj,1)
                  zgdrho = MAX( 0.e0, zgdrho )                               ! only if shelf is denser than deep
                  !
                  !                                                          ! bbl transport (down-slope direction)
                  utr_bbl(ji,jj) = e2u(ji,jj) * e3u_bbl_0(ji,jj) * zgbbl * zgdrho * REAL( mgrhu(ji,jj) )
                  !
                  !                                                  ! j-direction
                  !  down-slope T-point j/k-index (deep)  &   of the up  -slope T-point j/k-index (shelf)
                  ijd  = jj + MAX( 0, mgrhv(ji,jj) )
                  ijs  = jj + 1 - MAX( 0, mgrhv(ji,jj) )
                  !
                  ikvd = mbkv_d(ji,jj)
                  ikvs = mbkv(ji,jj)
                  !
                  za = zab(ji,jj+1,jp_tem) + zab(ji,jj,jp_tem)               ! 2*(alpha,beta) at v-point
                  zb = zab(ji,jj+1,jp_sal) + zab(ji,jj,jp_sal)
                  !                                                          !   masked bottom density gradient
                  zgdrho = 0.5 * (  za * ( zts(ji,ijd,jp_tem) - zts(ji,ijs,jp_tem) )    &
                     &            - zb * ( zts(ji,ijd,jp_sal) - zts(ji,ijs,jp_sal) )  ) * vmask(ji,jj,1)
                  zgdrho = MAX( 0.e0, zgdrho )                               ! only if shelf is denser than deep
                  !
                  !                                                          ! bbl transport (down-slope direction)
                  vtr_bbl(ji,jj) = e1v(ji,jj) * e3v_bbl_0(ji,jj) * zgbbl * zgdrho * REAL( mgrhv(ji,jj) )
               END DO
            END DO
         END SELECT
         !
      ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop( 'bbl')
      !
   END SUBROUTINE bbl


   SUBROUTINE tra_bbl_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_bbl_init  ***
      !!
      !! ** Purpose :   Initialization for the bottom boundary layer scheme.
      !!
      !! ** Method  :   Read the nambbl namelist and check the parameters
      !!              called by nemo_init at the first timestep (kit000)
      !!----------------------------------------------------------------------
      INTEGER ::   ji, jj               ! dummy loop indices
      INTEGER ::   ii0, ii1, ij0, ij1   ! local integer
      INTEGER ::   ios                  !   -      -
      REAL(wp), POINTER, DIMENSION(:,:) :: zmbk
      !!
      NAMELIST/nambbl/ nn_bbl_ldf, nn_bbl_adv, rn_ahtbbl, rn_gambbl, ln_kriteria
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start( 'tra_bbl_init')
      !
      CALL wrk_alloc( jpi, jpj, zmbk )
      !

      REWIND( numnam_ref )              ! Namelist nambbl in reference namelist : Bottom boundary layer scheme
      READ  ( numnam_ref, nambbl, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nambbl in reference namelist', lwp )

      REWIND( numnam_cfg )              ! Namelist nambbl in configuration namelist : Bottom boundary layer scheme
      READ  ( numnam_cfg, nambbl, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nambbl in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, nambbl )
      !
      l_bbl = .TRUE.                 !* flag to compute bbl coef and transport
      !
      IF(lwp) THEN                   !* Parameter control and print
         WRITE(numout,*)
         WRITE(numout,*) 'tra_bbl_init : bottom boundary layer initialisation'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '       Namelist nambbl : set bbl parameters'
         WRITE(numout,*) '          diffusive bbl (=1)   or not (=0)    nn_bbl_ldf = ', nn_bbl_ldf
         WRITE(numout,*) '          advective bbl (=1/2) or not (=0)    nn_bbl_adv = ', nn_bbl_adv
         WRITE(numout,*) '          diffusive bbl coefficient           rn_ahtbbl  = ', rn_ahtbbl, ' m2/s'
         WRITE(numout,*) '          advective bbl coefficient           rn_gambbl  = ', rn_gambbl, ' s'
         IF( nn_bbl_adv == 1 )    WRITE(numout,*) '       * Advective BBL using upper velocity'
         IF( nn_bbl_adv == 2 )    WRITE(numout,*) '       * Advective BBL using velocity = F( delta rho)'
      ENDIF

      !                              ! allocate trabbl arrays
      IF( tra_bbl_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'tra_bbl_init : unable to allocate arrays' )
     
!JMM no more in 3.6 ???
      IF( nn_eos /= 0 )   CALL ctl_stop ( ' bbl parameterisation requires eos = 0. We stop.' )

      !                             !* vertical index of  "deep" bottom u- and v-points
      DO jj = 1, jpjm1                    ! (the "shelf" bottom k-indices are mbku and mbkv)
         DO ji = 1, jpim1
            mbku_d(ji,jj) = MAX(  mbkt(ji+1,jj  ) , mbkt(ji,jj)  )   ! >= 1 as mbkt=1 over land
            mbkv_d(ji,jj) = MAX(  mbkt(ji  ,jj+1) , mbkt(ji,jj)  )
         END DO
      END DO
      ! convert into REAL to use lbc_lnk ; impose a min value of 1 as a zero can be set in lbclnk
      zmbk(:,:) = REAL( mbku_d(:,:), wp )   ;   CALL lbc_lnk(zmbk,'U',1.)   ;   mbku_d(:,:) = MAX( INT( zmbk(:,:) ), 1 )
      zmbk(:,:) = REAL( mbkv_d(:,:), wp )   ;   CALL lbc_lnk(zmbk,'V',1.)   ;   mbkv_d(:,:) = MAX( INT( zmbk(:,:) ), 1 )

                                        !* sign of grad(H) at u- and v-points
      mgrhu(jpi,:) = 0   ;   mgrhu(:,jpj) = 0   ;   mgrhv(jpi,:) = 0   ;   mgrhv(:,jpj) = 0
      DO jj = 1, jpjm1
         DO ji = 1, jpim1
! { DRAKKAR : ln_kriteria
            IF ( ln_kriteria ) THEN
              mgrhu(ji,jj) = mbkt(ji+1,jj) - mbkt(ji,jj)
              mgrhv(ji,jj) = mbkt(ji,jj+1) - mbkt(ji,jj)
              IF ( mgrhu(ji,jj) /= 0 ) THEN 
                  mgrhu(ji,jj) = INT(  SIGN( 1.e0, FLOAT(mbkt(ji+1,jj) - mbkt(ji,jj)) )) 
              ENDIF
              IF ( mgrhv(ji,jj) /= 0 ) THEN 
                  mgrhv(ji,jj) = INT(  SIGN( 1.e0, FLOAT(mbkt(ji,jj+1) - mbkt(ji,jj)) ))  
              ENDIF
            ELSE
            mgrhu(ji,jj) = INT(  SIGN( 1.e0, gdept_0(ji+1,jj,mbkt(ji+1,jj)) - gdept_0(ji,jj,mbkt(ji,jj)) )  )
            mgrhv(ji,jj) = INT(  SIGN( 1.e0, gdept_0(ji,jj+1,mbkt(ji,jj+1)) - gdept_0(ji,jj,mbkt(ji,jj)) )  )
            ENDIF
         END DO
      END DO

      DO jj = 1, jpjm1              !* bbl thickness at u- (v-) point
         DO ji = 1, jpim1                 ! minimum of top & bottom e3u_0 (e3v_0)
            e3u_bbl_0(ji,jj) = MIN( e3u_0(ji,jj,mbkt(ji+1,jj  )), e3u_0(ji,jj,mbkt(ji,jj)) )
            e3v_bbl_0(ji,jj) = MIN( e3v_0(ji,jj,mbkt(ji  ,jj+1)), e3v_0(ji,jj,mbkt(ji,jj)) )
         END DO
      END DO
      CALL lbc_lnk( e3u_bbl_0, 'U', 1. )   ;   CALL lbc_lnk( e3v_bbl_0, 'V', 1. )      ! lateral boundary conditions

      !                             !* masked diffusive flux coefficients
      ahu_bbl_0(:,:) = rn_ahtbbl * e2u(:,:) * e3u_bbl_0(:,:) / e1u(:,:)  * umask(:,:,1)
      ahv_bbl_0(:,:) = rn_ahtbbl * e1v(:,:) * e3v_bbl_0(:,:) / e2v(:,:)  * vmask(:,:,1)


      IF( cp_cfg == "orca" ) THEN   !* ORCA configuration : regional enhancement of ah_bbl
         !
         SELECT CASE ( jp_cfg )
         CASE ( 2 )                          ! ORCA_R2
            ij0 = 102   ;   ij1 = 102              ! Gibraltar enhancement of BBL
            ii0 = 139   ;   ii1 = 140
            ahu_bbl_0(mi0(ii0):mi1(ii1),mj0(ij0):mj1(ij1)) =  4.e0*ahu_bbl_0(mi0(ii0):mi1(ii1),mj0(ij0):mj1(ij1))
            ahv_bbl_0(mi0(ii0):mi1(ii1),mj0(ij0):mj1(ij1)) =  4.e0*ahv_bbl_0(mi0(ii0):mi1(ii1),mj0(ij0):mj1(ij1))
            !
            ij0 =  88   ;   ij1 =  88              ! Red Sea enhancement of BBL
            ii0 = 161   ;   ii1 = 162
            ahu_bbl_0(mi0(ii0):mi1(ii1),mj0(ij0):mj1(ij1)) = 10.e0*ahu_bbl_0(mi0(ii0):mi1(ii1),mj0(ij0):mj1(ij1))
            ahv_bbl_0(mi0(ii0):mi1(ii1),mj0(ij0):mj1(ij1)) = 10.e0*ahv_bbl_0(mi0(ii0):mi1(ii1),mj0(ij0):mj1(ij1))
            !
         CASE ( 4 )                          ! ORCA_R4
            ij0 =  52   ;   ij1 =  52              ! Gibraltar enhancement of BBL
            ii0 =  70   ;   ii1 =  71
            ahu_bbl_0(mi0(ii0):mi1(ii1),mj0(ij0):mj1(ij1)) =  4.e0*ahu_bbl_0(mi0(ii0):mi1(ii1),mj0(ij0):mj1(ij1))
            ahv_bbl_0(mi0(ii0):mi1(ii1),mj0(ij0):mj1(ij1)) =  4.e0*ahv_bbl_0(mi0(ii0):mi1(ii1),mj0(ij0):mj1(ij1))
         END SELECT
         !
      ENDIF
      !
      CALL wrk_dealloc( jpi, jpj, zmbk )
      !
      IF( nn_timing == 1 )  CALL timing_stop( 'tra_bbl_init')
      !
   END SUBROUTINE tra_bbl_init

#else
   !!----------------------------------------------------------------------
   !!   Dummy module :                      No bottom boundary layer scheme
   !!----------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_trabbl = .FALSE.   !: bbl flag
CONTAINS
   SUBROUTINE tra_bbl_init               ! Dummy routine
   END SUBROUTINE tra_bbl_init
   SUBROUTINE tra_bbl( kt )              ! Dummy routine
      WRITE(*,*) 'tra_bbl: You should not have seen this print! error?', kt
   END SUBROUTINE tra_bbl
#endif

   !!======================================================================
END MODULE trabbl
