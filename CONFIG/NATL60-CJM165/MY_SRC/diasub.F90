MODULE diasub
   !!======================================================================
   !!                       ***  MODULE  diasub  ***
   !! Ocean diagnostics: subduction rates
   !!======================================================================
   !! History :  OPA  !  2006-06  (P. Karleskind)  Original code
   !!   NEMO     3.2  !  2011-09  (P. Karleskind)  Nemo version
   !!----------------------------------------------------------------------
#if defined key_diasub
   !!----------------------------------------------------------------------
   !!   'key_diasub' :                              subduction diag.
   !!----------------------------------------------------------------------
   !!   
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! MPP library
   USE iom             ! I/O library
   USE zdfmxl, ONLY :  nmln  !: number of level in the mixed layer
   USE zdfmxl, ONLY :  hmlp      !: MLD rho criterion

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dia_sub       ! routine called by step.F90
   PUBLIC   dia_sub_alloc ! routine called by nemogcm.F90

   LOGICAL , PUBLIC, PARAMETER ::   lk_diasub = .TRUE.  !: subduction flag

   ! note: following variables should move to local variables once iom_put is always used 
   INTEGER, PUBLIC, PARAMETER ::   jpsub_xad     =  1   !: x- horizontal advection
   INTEGER, PUBLIC, PARAMETER ::   jpsub_yad     =  2   !: y- horizontal advection
   INTEGER, PUBLIC, PARAMETER ::   jpsub_zad     =  3   !: z- vertical   advection
   INTEGER, PUBLIC, PARAMETER ::   jpsub_mld     =  4   !: mixed layer depth
   INTEGER, PUBLIC, PARAMETER ::   jpsub_xlf     =  5   !: x- lateral       diffusion
   INTEGER, PUBLIC, PARAMETER ::   jpsub_ylf     =  6   !: y- lateral       diffusion
   INTEGER, PUBLIC, PARAMETER ::   jpsub_zlf     =  7  !: y- lateral       diffusion
   INTEGER, PUBLIC, PARAMETER ::   jpsub_zdf     =  8  !: vertical diffusion (Kz)
#if defined key_trcldf_eiv
   INTEGER, PUBLIC, PARAMETER ::   jpsub_xei     =  9   !: x- horiz. EIV advection
   INTEGER, PUBLIC, PARAMETER ::   jpsub_yei     =  10   !: y- horiz. EIV advection
   INTEGER, PUBLIC, PARAMETER ::   jpsub_zei     =  11   !: z- vert.  EIV advection
   INTEGER, PUBLIC, PARAMETER ::   jptrsub       =  11  !: number of subduction diagnostics type
#else
   INTEGER, PUBLIC, PARAMETER ::   jptrsub       =  8  !: number of subduction diagnostics type
#endif
   !
   INTEGER , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::  nkind      !: index of subduction for w
   INTEGER , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::  nindui   !: index of subduction for u
   INTEGER , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::  nindvj   !: index of subduction for v
   INTEGER , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::  nmlb    !:before mld index
   INTEGER , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::  minku, maxku
   INTEGER , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::  minkv, maxkv
   INTEGER , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::  minh, maxh
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::  tsub     !: subduction rate
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::  hmlpb     !: before mld
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::  gmlh     !: mld term
   !
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::  vht      !: vertical velocity subduction 
   !

   !! * Substitutions
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO Consortium (2011)
   !! $Id: diasub.F90 2715 2011-03-30 15:58:35Z rblod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION dia_sub_alloc()
      !!---------------------------------------------------------------------
      !
      ALLOCATE( nkind(jpi,jpj), nindvj(jpi,jpj,jpk), nindui(jpi,jpj,jpk), nmlb(jpi,jpj),   &
            &   minku(jpi,jpj), maxku(jpi,jpj), minkv(jpi,jpj), maxkv(jpi,jpj),            & 
            &   minh(jpi,jpj), maxh(jpi,jpj), hmlpb(jpi,jpj), gmlh(jpi,jpj,jpk),            &
            &   tsub(jpi,jpj,4),                                                          &
            &                                                         STAT=dia_sub_alloc)
      !
      IF( lk_mpp           )   CALL mpp_sum ( dia_sub_alloc )
      IF(dia_sub_alloc /= 0)   CALL ctl_warn('dia_sub_alloc: failed to allocate arrays.')
      !
   END FUNCTION dia_sub_alloc


   SUBROUTINE dia_sub( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dia_sub  ***
      !!
      !! ** Purpose : Computes
      !!
      !! ** Method : 
      !!-------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
      !!
      INTEGER  ::  ji, jj, jk          ! dummy loop arguments
      INTEGER  ::  ierr, ikn, ikb, ikni, iknj
      REAL(wp) ::  zau, zav

      !!----------------------------------------------------------------------

      IF( kt == nit000 ) THEN

         IF (lwp ) WRITE(numout,*)' DIASUB routine : initialize subduction computation'
         !                                      ! allocate dia_sub array
         IF( dia_sub_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'diasub : unable to allocate standard arrays' )
         !
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dia_sub : diagnostics of the subduction rates'
         IF(lwp) WRITE(numout,*) '~~~~~~~ '
         IF(lwp) WRITE(numout,*)
         gmlh(:,:,:) = 0._wp
         tsub(:,:,:) = 0._wp 
      ENDIF

      DO jj = 2, jpj
        DO ji = 2, jpi
           ikn = nmln(ji,jj)  ;   ikni =  nmln(ji-1,jj)   ;    iknj = nmln(ji,jj-1) 
           !                                                   !   w-term of subduction rate
           tsub(ji,jj,jpsub_zad) = wn(ji,jj,ikn) * e1t(ji,jj) * e2t(ji,jj)
           !
           IF( wn(ji,jj,ikn) < 0._wp ) THEN   ;  nkind(ji,jj) = ikn - 1
           ELSE                               ;  nkind(ji,jj) = ikn
           ENDIF
           !                                                   !   u-advection term
           minku(ji,jj) = min( ikn, ikni )
           maxku(ji,jj) = max( ikn, ikni ) - 1
           DO jk = minku(ji,jj),maxku(ji,jj)
              zau = umask(ji-1,jj,jk)      &
               &    * sign( 1, ikni - ikn )  &
               &    * un(ji-1,jj,jk) * fse3t(ji-1,jj,jk) * e2u(ji-1,jj)                      
              tsub(ji,jj,jpsub_xad) =  tsub(ji,jj,jpsub_xad) + zau 
              !
              IF( un(ji-1,jj,jk) < 0._wp ) THEN  ;   nindui(ji,jj,jk) = ji
              ELSE                               ;   nindui(ji,jj,jk) = ji - 1
              ENDIF
           END DO
           !                                                   !  v-advection term
           minkv(ji,jj) = min( ikn , iknj )
           maxkv(ji,jj) = max( ikn , iknj ) - 1
           DO jk = minkv(ji,jj), maxkv(ji,jj)
              zav = vmask(ji,jj-1,jk)      &
               &     * sign( 1, iknj - ikn )  &
               &     * vn(ji,jj-1,jk) * fse3t(ji,jj-1,jk) * e1v(ji,jj-1)
              tsub(ji,jj,jpsub_yad) =  tsub(ji,jj,jpsub_yad) + zav 
              !
              IF( vn(ji,jj-1,jk) < 0._wp ) THEN  ;   nindvj(ji,jj,jk) = jj
              ELSE                               ;   nindvj(ji,jj,jk) = jj - 1
              ENDIF
            END DO
        END DO
      END DO
      !                                                   !  mld term rate
      IF( kt > nit000 ) THEN
        DO jj = 2, jpj   
          DO ji = 2, jpi  
            ikn = nmln(ji,jj)   ;     ikb = nmlb(ji,jj) 
            IF( ikn /= ikb ) THEN
               minh(ji,jj) = min( ikn, ikb ) 
               maxh(ji,jj) = max( ikn, ikb ) - 1 
               !                                           !  ml variablity term 
               gmlh(ji,jj,: ) = 0._wp
               DO jk = minh(ji,jj), maxh(ji,jj)
                  gmlh(ji,jj,jk) =  tmask(ji,jj,jk)       &
                   &                 * sign( 1, ikn - ikb )  &
                   &                 * e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,jk) / rdt
                  tsub(ji,jj,jpsub_mld) = tsub(ji,jj,jpsub_mld) - gmlh(ji,jj,jk)
               END DO
            ENDIF
          END DO
        END DO
      ENDIF
      !
      CALL iom_put( "subrate_xad" , tsub(:,:,jpsub_xad) )    !  U-term
      CALL iom_put( "subrate_yad" , tsub(:,:,jpsub_yad) )    !  V-term
      CALL iom_put( "subrate_zad" , tsub(:,:,jpsub_zad) )    !  W-term
      CALL iom_put( "subrate_mld" , tsub(:,:,jpsub_mld) )    !  MLD-term

   END SUBROUTINE dia_sub

#else
   !!----------------------------------------------------------------------
   !!   Default option :                                       Empty module
   !!----------------------------------------------------------------------
   LOGICAL , PUBLIC, PARAMETER ::   lk_diasub = .FALSE.  !: thermocline-20d depths flag
CONTAINS
   SUBROUTINE dia_sub( kt )         ! Empty routine
      !WRITE(*,*) 'dia_sub: You should not have seen this print! error?', kt
   END SUBROUTINE dia_sub
#endif

   !!======================================================================
END MODULE diasub
