MODULE shapiro 
   !!==============================================================================
   !!                       ***  MODULE shapiro   ***
   !! spatial filtering of input field
   !!==============================================================================
   !! History :       !  09-08  (S. Cailleau ) from N. Ferry
                      !  09-09  (C. Regnier  ) corrections
                      !  04-10  (J.M. Molines) module and nemo standard
   !!----------------------------------------------------------------------
   !! * Modules used
   USE in_out_manager
   USE dom_oce         ! ocean space and time domain
   USE wrk_nemo        ! work arrays
   USE timing          ! Timing
   USE lbclnk

   
   IMPLICIT NONE
   PRIVATE

   PUBLIC Shapiro_1D  ! use by sbcblk_core  and sbcssr

  CONTAINS

  SUBROUTINE Shapiro_1D(ptabin, kiter, cd_overlap, ptabout) !GIG
      !!----------------------------------------------------------------------
      !!                  ***  routine Shapiro_1D  ***
      !!
      !! ** Purpose :  Multiple application (kiter) of a shapiro filter
      !!               on ptabin to produce ptabout.
      !!
      !! ** Method  :   
      !!
      !! ** Action  :   ptabout filtered output from ptabin
      !!
      !!----------------------------------------------------------------------
      INTEGER,                  INTENT(IN)  :: kiter      ! number of iterations to perform
      REAL(wp), DIMENSION(:,:), INTENT(IN)  :: ptabin     ! input array
      CHARACTER(len=*),         INTENT(IN)  :: cd_overlap ! = one of MERCA_GLOB, REGULAR_GLOB, ORCA_GLOB (??)
      REAL(wp), DIMENSION(:,:), INTENT(OUT) :: ptabout    ! output array

      ! * Local variable
      INTEGER                               :: ji, jj, jn ! dummy loop index
      REAL(wp), POINTER, DIMENSION(:,:)     :: zvarout    ! working array
      REAL(wp), PARAMETER                   :: rp_aniso_diff_XY=2.25 !  anisotrope case (???)
                                                          ! Empirical value for 140 iterations
                                                          ! for an anisotropic ratio of 1.5.
                                                          ! (re ???)
      REAL(wp)                              :: zalphax    ! weight coeficient (x direction)
      REAL(wp)                              :: zalphay    ! weight coeficient (y direction)
      REAL(wp)                              :: znum       ! numerator
      REAL(wp)                              :: zden       ! denominator
!
!------------------------------------------------------------------------------
!
      IF( nn_timing == 1 )  CALL timing_start('Shapiro')
      !
      CALL wrk_alloc( jpi,jpj, zvarout )

!     Global ocean case
      IF (( cd_overlap == 'MERCA_GLOB' )   .OR.   &
          ( cd_overlap == 'REGULAR_GLOB' ) .OR.   &
          ( cd_overlap == 'ORCA_GLOB' )) THEN
       ptabout(:,:) = ptabin(:,:)
       zvarout(:,:) = ptabout(:,:)  ! ptabout intent out ???

       zalphax=1./2.
       zalphay=1./2.

!  Dx/Dy=rp_aniso_diff_XY  , D_ = vitesse de diffusion
!  140 passes du fitre, Lx/Ly=1.5, le rp_aniso_diff_XY correspondant est:

       IF ( rp_aniso_diff_XY >=  1. ) zalphay=zalphay/rp_aniso_diff_XY
       IF ( rp_aniso_diff_XY <   1. ) zalphax=zalphax*rp_aniso_diff_XY

        DO jn = 1,kiter
            DO jj = 2,jpjm1
               DO ji = 2,jpim1
                  ! We crop on the coast
                   znum = zvarout(ji,jj)   &
                          + 0.25*zalphax*(zvarout(ji-1,jj  )-zvarout(ji,jj))*tmask(ji-1,jj  ,1)  &
                          + 0.25*zalphax*(zvarout(ji+1,jj  )-zvarout(ji,jj))*tmask(ji+1,jj  ,1)  &
                          + 0.25*zalphay*(zvarout(ji  ,jj-1)-zvarout(ji,jj))*tmask(ji  ,jj-1,1)  &
                          + 0.25*zalphay*(zvarout(ji  ,jj+1)-zvarout(ji,jj))*tmask(ji  ,jj+1,1)
                   ptabout(ji,jj)=znum*tmask(ji,jj,1)+ptabin(ji,jj)*(1.-tmask(ji,jj,1))
                ENDDO  ! end loop jj
            ENDDO  ! end loop ji
!
!
!           Periodical condition in case of cd_overlap (global ocean)
!           - on a mercator projection grid we consider that singular point at poles
!             are a mean of the values at points of the previous latitude
!           - on ORCA and regular grid we copy the values at points of the previous latitude
            IF ( cd_overlap == 'MERCAT_GLOB' ) THEN
!GIG case unchecked  ! JMM for sure not valid in MPP (BUG)
               ptabout(1,1) = SUM(ptabout(:,2)) / jpi
               ptabout(jpi,jpj) = SUM(ptabout(:,jpj-1)) / jpi
            ELSE
               CALL lbc_lnk(ptabout, 'T', 1.) ! Boundary condition
            ENDIF
            zvarout(:,:) = ptabout(:,:)
         ENDDO  ! end loop jn
      ENDIF

      CALL wrk_dealloc( jpi,jpj, zvarout )
      IF( nn_timing == 1 )  CALL timing_stop('Shapiro')
!
END SUBROUTINE Shapiro_1D     

END MODULE shapiro
