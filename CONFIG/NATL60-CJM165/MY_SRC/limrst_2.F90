MODULE limrst_2
   !!======================================================================
   !!                     ***  MODULE  limrst_2  ***
   !! Ice restart :  write the ice restart file
   !!======================================================================
   !! History :  2.0  !  01-04  (C. Ethe, G. Madec)  Original code
   !!                 !  06-07  (S. Masson)  use IOM for restart read/write
   !!            3.3  !  09-05  (G.Garric) addition of the lim2_evp case
   !!----------------------------------------------------------------------
#if defined key_lim2
   !!----------------------------------------------------------------------
   !!   'key_lim2' :                                  LIM 2.0 sea-ice model
   !!----------------------------------------------------------------------
   !!   lim_rst_opn_2   : open ice restart file
   !!   lim_rst_write_2 : write of the ice restart file 
   !!   lim_rst_read_2  : read  the ice restart file 
   !!----------------------------------------------------------------------
   USE dom_oce          ! ocean domain
   USE ice_2            ! LIM-2: sea-ice variables
   USE sbc_oce          ! Surface Boundary Condition: ocean
   USE sbc_ice          ! Surface Boundary Condition: sea-ice
   USE in_out_manager   ! I/O manager
   USE iom              ! I/O library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   lim_rst_opn_2     ! routine called by sbcice_lim_2.F90
   PUBLIC   lim_rst_write_2   ! routine called by sbcice_lim_2.F90
   PUBLIC   lim_rst_read_2    ! routine called by iceini_2.F90

   LOGICAL, PUBLIC ::   lrst_ice         !: logical to control the ice restart write 
   INTEGER, PUBLIC ::   numrir, numriw   !: logical unit for ice restart (read and write)

   !!----------------------------------------------------------------------
   !! NEMO/LIM2 3.3 , UCL - NEMO Consortium (2010)
   !! $Id: limrst_2.F90 5341 2015-06-03 14:59:46Z davestorkey $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE lim_rst_opn_2( kt )
      !!----------------------------------------------------------------------
      !!                    ***  lim_rst_opn_2  ***
      !!
      !! ** purpose  :   output of sea-ice variable in a netcdf file
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt       ! number of iteration
      !
      CHARACTER(LEN=20)   ::   clkt     ! ocean time-step deine as a character
      CHARACTER(LEN=50)   ::   clname   ! ice output restart file name
      CHARACTER(len=150)  ::   clpath   ! full path to ice output restart file
      !!----------------------------------------------------------------------
      !
      IF( kt == nit000 )   lrst_ice = .FALSE.   ! default definition
      
      ! to get better performances with NetCDF format:
      ! we open and define the ice restart file one ice time step before writing the data (-> at nitrst - 2*nn_fsbc + 1)
      ! except if we write ice restart files every ice time step or if an ice restart file was writen at nitend - 2*nn_fsbc + 1
      IF( kt == nitrst - 2*nn_fsbc + 1 .OR. nstock == nn_fsbc .OR. ( kt == nitend - nn_fsbc + 1 .AND. .NOT. lrst_ice ) ) THEN
         IF( nitrst <= nitend .AND. nitrst > 0 ) THEN
            ! beware of the format used to write kt (default is i8.8, that should be large enough...)
            IF( nitrst > 99999999 ) THEN   ;   WRITE(clkt, *       ) nitrst
            ELSE                           ;   WRITE(clkt, '(i8.8)') nitrst
            ENDIF
            ! create the file
!{ DRAKKAR : use simpler name for restart files ( defined in iceini_2.F90 )
            !clname = TRIM(cexper)//"_"//TRIM(ADJUSTL(clkt))//"_"//TRIM(cn_icerst_out)
            clname = TRIM(cn_icerst_out)
!}
            clpath = TRIM(cn_icerst_outdir) 
            IF( clpath(LEN_TRIM(clpath):) /= '/' ) clpath = TRIM(clpath)//'/' 
            IF(lwp) THEN
               WRITE(numout,*)
               SELECT CASE ( jprstlib )
               CASE ( jprstdimg )
                  WRITE(numout,*) '             open ice restart binary file: ',TRIM(clpath)//clname
               CASE DEFAULT
                  WRITE(numout,*) '             open ice restart NetCDF file: ',TRIM(clpath)//clname
               END SELECT
               IF( kt == nitrst - 2*nn_fsbc + 1 ) THEN   
                  WRITE(numout,*)         '             kt = nitrst - 2*nn_fsbc + 1 = ', kt,' date= ', ndastp
               ELSE   ;   WRITE(numout,*) '             kt = '                         , kt,' date= ', ndastp
               ENDIF
            ENDIF

            CALL iom_open( TRIM(clpath)//TRIM(clname), numriw, ldwrt = .TRUE., kiolib = jprstlib )
            lrst_ice = .TRUE.
         ENDIF
      ENDIF
      !
   END SUBROUTINE lim_rst_opn_2


   SUBROUTINE lim_rst_write_2( kt )
      !!----------------------------------------------------------------------
      !!                    ***  lim_rst_write_2  ***
      !!
      !! ** purpose  :   output of sea-ice variable in a netcdf file
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! number of iteration
      !
      INTEGER ::   iter   ! kt + nn_fsbc -1
      !!----------------------------------------------------------------------

      iter = kt + nn_fsbc - 1   ! ice restarts are written at kt == nitrst - nn_fsbc + 1

      IF( iter == nitrst ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'lim_rst_write_2 : write ice restart file  kt =', kt
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~~'
      ENDIF

      ! Write in numriw (if iter == nitrst)
      ! ------------------ 
      !                                                                     ! calendar control
      CALL iom_rstput( iter, nitrst, numriw, 'kt_ice', REAL( iter, wp) ) 
      
      CALL iom_rstput( iter, nitrst, numriw, 'hicif'      , hicif (:,:)   )      ! prognostic variables 
      CALL iom_rstput( iter, nitrst, numriw, 'hsnif'      , hsnif (:,:)   )
      CALL iom_rstput( iter, nitrst, numriw, 'frld'       , frld  (:,:)   )
      CALL iom_rstput( iter, nitrst, numriw, 'sist'       , sist  (:,:)   )
      CALL iom_rstput( iter, nitrst, numriw, 'tbif1'      , tbif  (:,:,1) )
      CALL iom_rstput( iter, nitrst, numriw, 'tbif2'      , tbif  (:,:,2) )
      CALL iom_rstput( iter, nitrst, numriw, 'tbif3'      , tbif  (:,:,3) )
      CALL iom_rstput( iter, nitrst, numriw, 'u_ice'      , u_ice (:,:)   )
      CALL iom_rstput( iter, nitrst, numriw, 'v_ice'      , v_ice (:,:)   )
      CALL iom_rstput( iter, nitrst, numriw, 'qstoif'     , qstoif(:,:)   )
      CALL iom_rstput( iter, nitrst, numriw, 'fsbbq'      , fsbbq (:,:)   )
#if ! defined key_lim2_vp
      CALL iom_rstput( iter, nitrst, numriw, 'stress1_i'  , stress1_i (:,:) )    ! EVP rheology
      CALL iom_rstput( iter, nitrst, numriw, 'stress2_i'  , stress2_i (:,:) )
      CALL iom_rstput( iter, nitrst, numriw, 'stress12_i' , stress12_i(:,:) )
#endif
      CALL iom_rstput( iter, nitrst, numriw, 'sxice'      , sxice (:,:)   )
      CALL iom_rstput( iter, nitrst, numriw, 'syice'      , syice (:,:)   )
      CALL iom_rstput( iter, nitrst, numriw, 'sxxice'     , sxxice(:,:)   )
      CALL iom_rstput( iter, nitrst, numriw, 'syyice'     , syyice(:,:)   )
      CALL iom_rstput( iter, nitrst, numriw, 'sxyice'     , sxyice(:,:)   )
      CALL iom_rstput( iter, nitrst, numriw, 'sxsn'       , sxsn  (:,:)   )
      CALL iom_rstput( iter, nitrst, numriw, 'sysn'       , sysn  (:,:)   )
      CALL iom_rstput( iter, nitrst, numriw, 'sxxsn'      , sxxsn (:,:)   )
      CALL iom_rstput( iter, nitrst, numriw, 'syysn'      , syysn (:,:)   )
      CALL iom_rstput( iter, nitrst, numriw, 'sxysn'      , sxysn (:,:)   )
      CALL iom_rstput( iter, nitrst, numriw, 'sxa'        , sxa   (:,:)   )
      CALL iom_rstput( iter, nitrst, numriw, 'sya'        , sya   (:,:)   )
      CALL iom_rstput( iter, nitrst, numriw, 'sxxa'       , sxxa  (:,:)   )
      CALL iom_rstput( iter, nitrst, numriw, 'syya'       , syya  (:,:)   )
      CALL iom_rstput( iter, nitrst, numriw, 'sxya'       , sxya  (:,:)   )
      CALL iom_rstput( iter, nitrst, numriw, 'sxc0'       , sxc0  (:,:)   )
      CALL iom_rstput( iter, nitrst, numriw, 'syc0'       , syc0  (:,:)   )
      CALL iom_rstput( iter, nitrst, numriw, 'sxxc0'      , sxxc0 (:,:)   )
      CALL iom_rstput( iter, nitrst, numriw, 'syyc0'      , syyc0 (:,:)   )
      CALL iom_rstput( iter, nitrst, numriw, 'sxyc0'      , sxyc0 (:,:)   )
      CALL iom_rstput( iter, nitrst, numriw, 'sxc1'       , sxc1  (:,:)   )
      CALL iom_rstput( iter, nitrst, numriw, 'syc1'       , syc1  (:,:)   )
      CALL iom_rstput( iter, nitrst, numriw, 'sxxc1'      , sxxc1 (:,:)   )
      CALL iom_rstput( iter, nitrst, numriw, 'syyc1'      , syyc1 (:,:)   )
      CALL iom_rstput( iter, nitrst, numriw, 'sxyc1'      , sxyc1 (:,:)   )
      CALL iom_rstput( iter, nitrst, numriw, 'sxc2'       , sxc2  (:,:)   )
      CALL iom_rstput( iter, nitrst, numriw, 'syc2'       , syc2  (:,:)   )
      CALL iom_rstput( iter, nitrst, numriw, 'sxxc2'      , sxxc2 (:,:)   )
      CALL iom_rstput( iter, nitrst, numriw, 'syyc2'      , syyc2 (:,:)   )
      CALL iom_rstput( iter, nitrst, numriw, 'sxyc2'      , sxyc2 (:,:)   )
      CALL iom_rstput( iter, nitrst, numriw, 'sxst'       , sxst  (:,:)   )
      CALL iom_rstput( iter, nitrst, numriw, 'syst'       , syst  (:,:)   )
      CALL iom_rstput( iter, nitrst, numriw, 'sxxst'      , sxxst (:,:)   )
      CALL iom_rstput( iter, nitrst, numriw, 'syyst'      , syyst (:,:)   )
      CALL iom_rstput( iter, nitrst, numriw, 'sxyst'      , sxyst (:,:)   )
      
      IF( iter == nitrst ) THEN
         CALL iom_close( numriw )                         ! close the restart file
         lrst_ice = .FALSE.
      ENDIF
      !
   END SUBROUTINE lim_rst_write_2


   SUBROUTINE lim_rst_read_2
      !!----------------------------------------------------------------------
      !!                    ***  lim_rst_read_2  ***
      !!
      !! ** purpose  :   read of sea-ice variable restart in a netcdf file
      !!----------------------------------------------------------------------
      REAL(wp) ::   ziter
      INTEGER  ::   jlibalt = jprstlib
      LOGICAL  ::   llok
      INTEGER  ::   itest
      !!----------------------------------------------------------------------

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'lim_rst_read_2 : read ice NetCDF restart file'
         WRITE(numout,*) '~~~~~~~~~~~~~~'
      ENDIF

      IF ( jprstlib == jprstdimg ) THEN
        ! eventually read netcdf file (monobloc)  for restarting on different number of processors
        ! if {cn_icerst_in}.nc exists, then set jlibalt to jpnf90
        INQUIRE( FILE = TRIM(cn_icerst_indir)//'/'//TRIM(cn_icerst_in)//'.nc', EXIST = llok )
        IF ( llok ) THEN ; jlibalt = jpnf90  ; ELSE ; jlibalt = jprstlib ; ENDIF
      ENDIF

      CALL iom_open ( TRIM(cn_icerst_indir)//'/'//TRIM(cn_icerst_in), numrir, kiolib = jlibalt )

      CALL iom_get( numrir, 'kt_ice' , ziter )
      IF(lwp) WRITE(numout,*) '   read ice restart file at time step    : ', INT( ziter )
      IF(lwp) WRITE(numout,*) '   in any case we force it to nit000 - 1 : ', nit000 - 1

      !Control of date
      
      IF( ( nit000 - INT(ziter) ) /= 1 .AND. ABS( nrstdt ) == 1 )   &
         &     CALL ctl_stop( 'lim_rst_read ===>>>> : problem with nit000 in ice restart',  &
         &                   '   verify the file or rerun with the value 0 for the',        &
         &                   '   control of time parameter  nrstdt' )

      CALL iom_get( numrir, jpdom_autoglo, 'hicif' , hicif  )    
      CALL iom_get( numrir, jpdom_autoglo, 'hsnif' , hsnif  )    
      CALL iom_get( numrir, jpdom_autoglo, 'frld'  , frld   )    
      CALL iom_get( numrir, jpdom_autoglo, 'sist'  , sist   )    
      CALL iom_get( numrir, jpdom_autoglo, 'tbif1' , tbif(:,:,1) )    
      CALL iom_get( numrir, jpdom_autoglo, 'tbif2' , tbif(:,:,2) )    
      CALL iom_get( numrir, jpdom_autoglo, 'tbif3' , tbif(:,:,3) )    

      itest = iom_varid( numrir, 'u_ice', ldstop = .FALSE. )   ! test if the variable u_ice is included in the file
      IF( itest > 0 ) THEN   ! yes -> new restart files (from NEMO 3.2)
         CALL iom_get( numrir, jpdom_autoglo, 'u_ice', u_ice )  
         CALL iom_get( numrir, jpdom_autoglo, 'v_ice', v_ice )    
      ELSE                   ! no  -> old restart file with variable called [uv]i_ice (inroduced in NEMO 3.0)
         CALL iom_get( numrir, jpdom_autoglo, 'ui_ice', u_ice )  
         CALL iom_get( numrir, jpdom_autoglo, 'vi_ice', v_ice )    
      ENDIF

      CALL iom_get( numrir, jpdom_autoglo, 'qstoif'     , qstoif )    
      CALL iom_get( numrir, jpdom_autoglo, 'fsbbq'      , fsbbq  )    
#if ! defined key_lim2_vp
      CALL iom_get( numrir, jpdom_autoglo, 'stress1_i'  , stress1_i  )
      CALL iom_get( numrir, jpdom_autoglo, 'stress2_i'  , stress2_i  )
      CALL iom_get( numrir, jpdom_autoglo, 'stress12_i' , stress12_i )
#endif
      CALL iom_get( numrir, jpdom_autoglo, 'sxice'      , sxice  )
      CALL iom_get( numrir, jpdom_autoglo, 'syice'      , syice  )
      CALL iom_get( numrir, jpdom_autoglo, 'sxxice'     , sxxice )
      CALL iom_get( numrir, jpdom_autoglo, 'syyice'     , syyice )
      CALL iom_get( numrir, jpdom_autoglo, 'sxyice'     , sxyice )
      CALL iom_get( numrir, jpdom_autoglo, 'sxsn'       , sxsn   )
      CALL iom_get( numrir, jpdom_autoglo, 'sysn'       , sysn   )
      CALL iom_get( numrir, jpdom_autoglo, 'sxxsn'      , sxxsn  )
      CALL iom_get( numrir, jpdom_autoglo, 'syysn'      , syysn  )
      CALL iom_get( numrir, jpdom_autoglo, 'sxysn'      , sxysn  )
      CALL iom_get( numrir, jpdom_autoglo, 'sxa'        , sxa    )
      CALL iom_get( numrir, jpdom_autoglo, 'sya'        , sya    )
      CALL iom_get( numrir, jpdom_autoglo, 'sxxa'       , sxxa   )
      CALL iom_get( numrir, jpdom_autoglo, 'syya'       , syya   )
      CALL iom_get( numrir, jpdom_autoglo, 'sxya'       , sxya   )
      CALL iom_get( numrir, jpdom_autoglo, 'sxc0'       , sxc0   )
      CALL iom_get( numrir, jpdom_autoglo, 'syc0'       , syc0   )
      CALL iom_get( numrir, jpdom_autoglo, 'sxxc0'      , sxxc0  )
      CALL iom_get( numrir, jpdom_autoglo, 'syyc0'      , syyc0  )
      CALL iom_get( numrir, jpdom_autoglo, 'sxyc0'      , sxyc0  )
      CALL iom_get( numrir, jpdom_autoglo, 'sxc1'       , sxc1   )
      CALL iom_get( numrir, jpdom_autoglo, 'syc1'       , syc1   )
      CALL iom_get( numrir, jpdom_autoglo, 'sxxc1'      , sxxc1  )
      CALL iom_get( numrir, jpdom_autoglo, 'syyc1'      , syyc1  )
      CALL iom_get( numrir, jpdom_autoglo, 'sxyc1'      , sxyc1  )
      CALL iom_get( numrir, jpdom_autoglo, 'sxc2'       , sxc2   )
      CALL iom_get( numrir, jpdom_autoglo, 'syc2'       , syc2   )
      CALL iom_get( numrir, jpdom_autoglo, 'sxxc2'      , sxxc2  )
      CALL iom_get( numrir, jpdom_autoglo, 'syyc2'      , syyc2  )
      CALL iom_get( numrir, jpdom_autoglo, 'sxyc2'      , sxyc2  )
      CALL iom_get( numrir, jpdom_autoglo, 'sxst'       , sxst   )
      CALL iom_get( numrir, jpdom_autoglo, 'syst'       , syst   )
      CALL iom_get( numrir, jpdom_autoglo, 'sxxst'      , sxxst  )
      CALL iom_get( numrir, jpdom_autoglo, 'syyst'      , syyst  )
      CALL iom_get( numrir, jpdom_autoglo, 'sxyst'      , sxyst  )
      
      CALL iom_close( numrir )
      !
   END SUBROUTINE lim_rst_read_2

#else
   !!----------------------------------------------------------------------
   !!   Default option :       Empty module        NO LIM 2.0 sea-ice model
   !!----------------------------------------------------------------------
#endif

   !!======================================================================
END MODULE limrst_2
