MODULE trcwri
   !!======================================================================
   !!                       *** MODULE trcwri ***
   !!    TOP :   Output of passive tracers
   !!======================================================================
   !! History :   1.0  !  2009-05 (C. Ethe)  Original code
   !!----------------------------------------------------------------------
#if defined key_top && defined key_iomput
   !!----------------------------------------------------------------------
   !!   'key_top'                                           TOP models
   !!----------------------------------------------------------------------
   !! trc_wri_trc   :  outputs of concentration fields
   !!----------------------------------------------------------------------
   USE dom_oce     ! ocean space and time domain variables
   USE oce_trc     ! shared variables between ocean and passive tracers
   USE trc         ! passive tracers common variables 
   USE iom         ! I/O manager
   USE dianam      ! Output file name
   USE trcwri_pisces
   USE trcwri_cfc
   USE trcwri_c14b
   USE trcwri_my_trc
   USE diasub
   USE trcdiasub


   IMPLICIT NONE
   PRIVATE

   PUBLIC trc_wri      

   !! * Substitutions
#  include "top_substitute.h90"

CONTAINS

   SUBROUTINE trc_wri( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_wri  ***
      !! 
      !! ** Purpose :   output passive tracers fields and dynamical trends
      !!---------------------------------------------------------------------
      INTEGER, INTENT( in )     :: kt
      !
      INTEGER                   :: jn
      CHARACTER (len=20)        :: cltra
      CHARACTER (len=40)        :: clhstnam
      INTEGER ::   inum = 11            ! temporary logical unit
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('trc_wri')
      !
      IF( lk_offline .AND. kt == nittrc000 .AND. lwp ) THEN    ! WRITE root name in date.file for use by postpro
         CALL dia_nam( clhstnam, nn_writetrc,' ' )
         CALL ctl_opn( inum, 'date.file', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, numout, lwp, narea )
         WRITE(inum,*) clhstnam
         CLOSE(inum)
      ENDIF
      ! write the tracer concentrations in the file
      ! ---------------------------------------
      IF( lk_pisces  )   CALL trc_wri_pisces     ! PISCES 
      IF( lk_cfc     )   CALL trc_wri_cfc        ! surface fluxes of CFC
      IF( lk_c14b    )   CALL trc_wri_c14b       ! surface fluxes of C14
      IF( lk_my_trc  )   CALL trc_wri_my_trc     ! MY_TRC  tracers
      !
      IF( lk_diasub  )   CALL trc_wri_sub        ! Subduction diagnostics
      !
      IF( nn_timing == 1 )  CALL timing_stop('trc_wri')
      !
   END SUBROUTINE trc_wri

# if defined key_diasub

   SUBROUTINE trc_wri_sub
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE trc_wri_sub  ***
      !!
      !! ** Purpose :   output of passive tracer : advection-diffusion  subduction subduction
      !!
      !!----------------------------------------------------------------------
      CHARACTER(len=8), DIMENSION(jptrsub)   :: cltra1
      CHARACTER(len=20)  :: cltra, cltra2
      INTEGER            :: jn, jl, ji, jj, ik
      REAL(wp), DIMENSION(jpi,jpj)         :: zsed
      REAL(wp), DIMENSION(jpi,jpj,jptrsub) :: ztrsubtpoc
      !!----------------------------------------------------------------------

      DO jl = 1, jptrsub
         IF( jl == jpsub_xad ) cltra1(jl) = TRIM("xad_sub_")   ! x advection for tracer
         IF( jl == jpsub_yad ) cltra1(jl) = TRIM("yad_sub_")   ! y advection for tracer
         IF( jl == jpsub_zad ) cltra1(jl) = TRIM("zad_sub_")   ! z advection for tracer
         IF( jl == jpsub_mld ) cltra1(jl) = TRIM("mld_sub_")   ! mld for tracer
         IF( jl == jpsub_xlf ) cltra1(jl) = TRIM("xlf_sub_")   ! x lateral diffusion for tracer
         IF( jl == jpsub_ylf ) cltra1(jl) = TRIM("ylf_sub_")   ! y lateral diffusion for tracer
         IF( jl == jpsub_zlf ) cltra1(jl) = TRIM("zlf_sub_")   ! z lateral diffusion for tracer
         IF( jl == jpsub_zdf ) cltra1(jl) = TRIM("zdf_sub_")   ! z vertical diffusion for tracer
#if defined key_trcldf_eiv
         IF( jl == jpsub_xei ) cltra1(jl) = TRIM("xei_sub_")   ! x gent velocity for tracer
         IF( jl == jpsub_yei ) cltra1(jl) = TRIM("yei_sub_")   ! y gent velocity for tracer
         IF( jl == jpsub_zei ) cltra1(jl) = TRIM("zei_sub_")   ! z gent velocity for tracer
#endif
     ENDDO
     !
# if defined key_pisces
     DO jl = 1, jptrsub
        ! write the trends
        DO jn = 1, jptra
           IF(     jn == jpdic .OR. jn == jptal .OR. jn == jpoxy .OR. jn == jpdoc .OR. jn == jpcal  &
            & .OR. jn == jpno3 .OR. jn == jppo4 .OR. jn == jpsil .OR. jn == jpfer .OR. jn == jpnh4 ) THEN
              cltra = TRIM(cltra1(jl))//TRIM(ctrcnm(jn))
              CALL iom_put( cltra, trsub(:,:,jn,jl) )
           ENDIF
        END DO
        !
        ! Sum of subduction of phy, phy2, zoo, zoo2, POC et GOC
        ztrsubtpoc(:,:,jl) = trsub(:,:,jpphy,jl) + trsub(:,:,jpdia,jl)  &
           &               + trsub(:,:,jpzoo,jl) + trsub(:,:,jpmes,jl)  &
           &               + trsub(:,:,jppoc,jl) + trsub(:,:,jpgoc,jl)
        cltra2 = "TOC"
        cltra  = TRIM(cltra1(jl))//TRIM(cltra2)
        CALL iom_put( cltra, ztrsubtpoc(:,:,jl))
     END DO

     ! Sedimentation a la base de la couche de mélane
     cltra = TRIM("sedsed")
     DO jj = 1, jpj
        DO ji = 1, jpi
           ik = nmln(ji,jj)
           zsed(ji,jj) = ( sinking(ji,jj,ik) + sinking2(ji,jj,ik) ) * tmask(ji,jj,ik)
        END DO
     END DO
     CALL iom_put( cltra, zsed(:,:) )
#else
     DO jl = 1, jptrsub
        ! write the trends
        DO jn = 1, jptra
             cltra = TRIM(cltra1(jl))//TRIM(ctrcnm(jn))
             CALL iom_put( cltra, trsub(:,:,jn,jl) )
        END DO
        !
     END DO
#endif
      !
      ! une fois ecrit, trsub est remis a 0
      DO jl = 1, jptrsub
         DO jn = 1, jptra
            trsub(:,:,jn,jl) = 0.
         END DO
      END DO
      !
   END SUBROUTINE trc_wri_sub

# endif

#else
   !!----------------------------------------------------------------------
   !!  Dummy module :                                     No passive tracer
   !!----------------------------------------------------------------------
   PUBLIC trc_wri
CONTAINS
   SUBROUTINE trc_wri( kt )                     ! Empty routine   
   INTEGER, INTENT(in) :: kt
   END SUBROUTINE trc_wri
#endif

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcwri.F90 3750 2013-01-14 16:25:10Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!======================================================================
END MODULE trcwri
