MODULE trcini_my_trc
   !!======================================================================
   !!                         ***  MODULE trcini_my_trc  ***
   !! TOP :   initialisation of the MY_TRC tracers
   !!======================================================================
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec) Original code
   !!----------------------------------------------------------------------
#if defined key_my_trc
   !!----------------------------------------------------------------------
   !!   'key_my_trc'                                               CFC tracers
   !!----------------------------------------------------------------------
   !! trc_ini_my_trc   : MY_TRC model initialisation
   !!----------------------------------------------------------------------
   USE par_trc         ! TOP parameters
   USE oce_trc
   USE trc
   USE trcsms_my_trc

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_ini_my_trc   ! called by trcini.F90 module

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcini_my_trc.F90 5385 2015-06-09 13:50:42Z cetlod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_ini_my_trc
      !!----------------------------------------------------------------------
      !!                     ***  trc_ini_my_trc  ***  
      !!
      !! ** Purpose :   initialization for MY_TRC model
      !!
      !!----------------------------------------------------------------------

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_ini_my_trc: passive tracer unit vector'
      IF(lwp) WRITE(numout,*) ' To check conservation : '
      IF(lwp) WRITE(numout,*) '   1 - No sea-ice model '
      IF(lwp) WRITE(numout,*) '   2 - No runoff ' 
      IF(lwp) WRITE(numout,*) '   3 - precipitation and evaporation equal to 1 : E=P=1 ' 
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'
      
      IF( .NOT. ln_rsttr ) trn(:,:,:,jp_myt0:jp_myt1) = 0.
      !
   END SUBROUTINE trc_ini_my_trc

#else
   !!----------------------------------------------------------------------
   !!   Dummy module                                        No MY_TRC model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_ini_my_trc             ! Empty routine
   END SUBROUTINE trc_ini_my_trc
#endif

   !!======================================================================
END MODULE trcini_my_trc
