MODULE iom_rstdimg
   !!=====================================================================
   !!                    ***  MODULE  iom_rstdimg ***
   !! Input/Output manager :  Library to read input rstdimg files
   !!====================================================================
   !! History :  9.0  ! 06 09  (S. Masson) Original code
   !!--------------------------------------------------------------------

   !!--------------------------------------------------------------------
   !!   iom_open       : open a file read only
   !!   iom_close      : close a file or all files opened by iom
   !!   iom_get        : read a field (interfaced to several routines)
   !!   iom_gettime    : read the time axis kvid in the file
   !!   iom_varid      : get the id of a variable in a file
   !!   iom_rstput     : write a field in a restart file (interfaced to several routines)
   !!--------------------------------------------------------------------
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! MPP library
   USE dom_oce         ! ocean space and time domain
   USE lbclnk          ! lateal boundary condition / mpp exchanges
   USE iom_def         ! iom variables definitions

   IMPLICIT NONE
   PRIVATE

   PUBLIC iom_rstdimg_open, iom_rstdimg_close, iom_rstdimg_get, iom_rstdimg_rstput

   INTERFACE iom_rstdimg_get
      MODULE PROCEDURE iom_rstdimg_g0d, iom_rstdimg_g123d
   END INTERFACE
   INTERFACE iom_rstdimg_rstput
      MODULE PROCEDURE iom_rstdimg_rp0d, iom_rstdimg_rp123d
   END INTERFACE

   INTEGER, PARAMETER ::   jpvnl          = 32   ! variable name length
      
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: iom_rstdimg.F90 2715 2011-03-30 15:58:35Z rblod $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE iom_rstdimg_open( cdname, kiomid, ldwrt, ldok, kdompar )
      !!---------------------------------------------------------------------
      !!                   ***  SUBROUTINE  iom_open  ***
      !!
      !! ** Purpose :  open an input file read only (return 0 if not found)
      !!---------------------------------------------------------------------
      CHARACTER(len=*)       , INTENT(inout)           ::   cdname      ! File name
      INTEGER                , INTENT(  out)           ::   kiomid      ! iom identifier of the opened file
      LOGICAL                , INTENT(in   )           ::   ldwrt       ! read or write the file?
      LOGICAL                , INTENT(in   )           ::   ldok        ! check the existence 
      INTEGER, DIMENSION(2,5), INTENT(in   ), OPTIONAL ::   kdompar     ! domain parameters: 

      CHARACTER(LEN=100)                      ::   clinfo                     ! info character
      CHARACTER(LEN=100)                      ::   cltmp                      ! temporary character
      CHARACTER(LEN=10 )                      ::   clstatus                   ! status of opened file (REPLACE or NEW)
      INTEGER                                 ::   jv                         ! loop counter
      INTEGER                                 ::   istop                      ! temporary storage of nstop
      INTEGER                                 ::   idrst                      ! logical unit of the restart file
      INTEGER                                 ::   iln                        ! lengths of character
      INTEGER                                 ::   irecl8                     ! record length
      INTEGER                                 ::   ios                        ! IO status
      INTEGER                                 ::   irhd                       ! record of the header infos (if negative, new format with 
                                                                              ! layout array written starting irhd+1
      INTEGER                                 ::   ivnum                      ! number of variables
      INTEGER                                 ::   ishft                      ! counter shift
      INTEGER                                 ::   inx, iny, inz              ! x,y,z dimension of the variable
      INTEGER                                 ::   in0d, in1d, in2d, in3d     ! number of 0/1/2/3D variables
      INTEGER                                 ::   ipni, ipnj, ipnij, iarea   ! domain decomposition 
      INTEGER                                 ::   iiglo, ijglo               ! domain global size 
      INTEGER                                 ::   jl                         ! loop variable
      LOGICAL                                 ::   llclobber                  ! local definition of ln_clobber
      CHARACTER(LEN=jpvnl), DIMENSION(jpmax_vars) ::   clna0d, clna1d, clna2d, clna3d     ! name of 0/1/2/3D variables
      REAL(wp),             DIMENSION(jpmax_vars) ::   zval0d, zval1d, zval2d, zval3d     ! value of 0d variables or record
      !                                                                                   ! position for 1/2/3D variables
      !---------------------------------------------------------------------
      clinfo = '                    iom_rstdimg_open ~~~  '
      istop = nstop      ! store the actual value of nstop
      ios = 0            ! default definition
      kiomid = 0         ! default definition
      llclobber = ldwrt .AND. ln_clobber
      ! get a free unit
      idrst = get_unit()  ! get a free logical unit for the restart file
!!$#if defined key_agrif 
!!$      idrst = Agrif_Get_Unit()      
!!$#endif 
      ! Open the file...
      ! =============
      IF( ldok .AND. .NOT. llclobber ) THEN      ! Open existing file...
         ! find the record length
         OPEN( idrst, FILE = TRIM(cdname), FORM = 'unformatted', ACCESS = 'direct'   &
            &       , RECL = 8, STATUS = 'old', ACTION = 'read', IOSTAT = ios, ERR = 987 )
         READ( idrst, REC = 1, IOSTAT = ios, ERR = 987 ) irecl8
         CLOSE( idrst )
         ! Open the file with the appropriate record length and parameters
         IF( ldwrt ) THEN  ! ... in readwrite mode
            IF(lwp) WRITE(numout,*) TRIM(clinfo)//' open existing file: '//TRIM(cdname)//' in READWRITE mode'
            OPEN( idrst, FILE = TRIM(cdname), FORM = 'unformatted', ACCESS = 'direct'   &
               &       , RECL = irecl8, STATUS = 'old', ACTION = 'readwrite', IOSTAT = ios, ERR = 987 )
         ELSE              ! ... in read mode
            IF(lwp) WRITE(numout,*) TRIM(clinfo)//' open existing file: '//TRIM(cdname)//' in READ mode'
            OPEN( idrst, FILE = TRIM(cdname), FORM = 'unformatted', ACCESS = 'direct'   &
               &       , RECL = irecl8, STATUS = 'old', ACTION = 'read'     , IOSTAT = ios, ERR = 987 )
         ENDIF
      ELSE                                       ! the file does not exist (or we overwrite it)
         iln = INDEX( cdname, '.dimg' )
         IF( ldwrt ) THEN  ! the file should be open in readwrite mode so we create it...
            irecl8= MAX( kdompar(1,1) * kdompar(2,1) * wp, (  15 ) * 4 , jpnij * wp )  ! in newformat, layout array are saved 
                                                                                       ! at the end 1 per record
            IF( jpnij > 1 ) THEN
               WRITE(cltmp,'(a,a,i4.4,a)') cdname(1:iln-1), '_', narea, '.dimg'
               cdname = TRIM(cltmp)
            ENDIF
            IF(lwp) WRITE(numout,*) TRIM(clinfo)//' create new file: '//TRIM(cdname)//' in READWRITE mode'
            
            IF( llclobber ) THEN   ;   clstatus = 'REPLACE' 
            ELSE                   ;   clstatus = 'NEW'
            ENDIF
            OPEN( idrst, FILE = TRIM(cdname), FORM = 'UNFORMATTED', ACCESS = 'DIRECT'   &
               &       , RECL = irecl8, STATUS = TRIM(clstatus), ACTION = 'readwrite', IOSTAT = ios, ERR = 987 )
         ELSE              ! the file should be open for read mode so it must exist...
            CALL ctl_stop( TRIM(clinfo), ' should be impossible case...' )
         ENDIF
      ENDIF
      ! Performs checks on the file
      ! =============
      IF( ldok .AND. .NOT. llclobber ) THEN      ! old file
         READ( idrst, REC = 1   , IOSTAT = ios, ERR = 987 )              &
              &   irecl8, inx, iny, inz, in0d, in1d, in2d, in3d, irhd,   &
              &   ipni, ipnj, ipnij, iarea, iiglo, ijglo
         irhd=ABS(irhd)  ! to take into account new format with irhd < 0
         READ( idrst, REC = irhd, IOSTAT = ios, ERR = 987 )                       &
            &   clna0d(1:in0d), zval0d(1:in0d), clna1d(1:in1d), zval1d(1:in1d),   &
            &   clna2d(1:in2d), zval2d(1:in2d), clna3d(1:in3d), zval3d(1:in3d)
         clinfo = TRIM(clinfo)//' file '//TRIM(cdname)
         IF( iiglo /= jpiglo       )   CALL ctl_stop( TRIM(clinfo), 'Mismatch in global domain size in i direction' )
         IF( ijglo /= jpjglo       )   CALL ctl_stop( TRIM(clinfo), 'Mismatch in global domain size in j direction' )
         IF( ldwrt ) THEN
            IF( inx   /= kdompar(1,1) )   CALL ctl_stop( TRIM(clinfo), 'Mismatch in local domain size in i direction' )
            IF( iny   /= kdompar(2,1) )   CALL ctl_stop( TRIM(clinfo), 'Mismatch in local domain size in j direction' )
         ENDIF
         IF( inz   /= jpk          )   CALL ctl_stop( TRIM(clinfo), 'Mismatch in domain size in k direction' )
         IF( ipni  /= jpni         )   CALL ctl_stop( TRIM(clinfo), 'Processor splitting changed along I' )
         IF( ipnj  /= jpnj         )   CALL ctl_stop( TRIM(clinfo), 'Processor splitting changed along J' )
         IF( ipnij /= jpnij        )   CALL ctl_stop( TRIM(clinfo), 'Total number of processors changed' )
         IF( iarea /= narea        )   CALL ctl_stop( TRIM(clinfo), 'Mismatch in area numbering ...' )
      ENDIF
      ! fill file informations
      ! =============
      IF( istop == nstop ) THEN   ! no error within this routine
!does not work with some compilers         kiomid = MINLOC(iom_file(:)%nfid, dim = 1)
         kiomid = 0
         DO jl = jpmax_files, 1, -1
            IF( iom_file(jl)%nfid == 0 )   kiomid = jl
         ENDDO
         iom_file(kiomid)%name    = TRIM(cdname)
         iom_file(kiomid)%nfid    = idrst
         iom_file(kiomid)%iolib   = jprstdimg
         iom_file(kiomid)%iduld   = -1
         IF( ldok ) THEN      ! old file
            ! read variables informations from the file header
            IF(  TRIM(clna0d(1)) == 'no0d' )   in0d = 0
            IF(  TRIM(clna1d(1)) == 'no1d' )   in1d = 0
            IF(  TRIM(clna2d(1)) == 'no2d' )   in2d = 0
            IF(  TRIM(clna3d(1)) == 'no3d' )   in3d = 0
            ivnum = in0d + in1d + in2d + in3d
            iom_file(kiomid)%nvars            = ivnum
            iom_file(kiomid)%irec             = 2 + in1d + in2d + inz * in3d
            iom_file(kiomid)%luld(   1:ivnum) = .FALSE.
            iom_file(kiomid)%scf(    1:ivnum) = 1.
            ! scalar variable
            DO jv = 1, in0d
               iom_file(kiomid)%cn_var(jv) = TRIM(clna0d(jv))
               iom_file(kiomid)%nvid(  jv) = 1
               iom_file(kiomid)%ndims( jv) = 0
               iom_file(kiomid)%ofs(   jv) = zval0d(jv)   ! warning: trick... we use ofs to store the value
            END DO
            ! 1d variable
            ishft = in0d
            DO jv = 1, in1d
               iom_file(kiomid)%cn_var(    ishft + jv) = TRIM(clna1d(jv))
               iom_file(kiomid)%nvid(      ishft + jv) = zval1d(jv)
               iom_file(kiomid)%ndims(     ishft + jv) = 1
               iom_file(kiomid)%dimsz(1  , ishft + jv) = jpk
               iom_file(kiomid)%ofs(       ishft + jv) = 0.
            END DO
            ! 2d variable
            ishft = in0d + in1d
            DO jv = 1, in2d
               iom_file(kiomid)%cn_var(    ishft + jv) = TRIM(clna2d(jv))
               iom_file(kiomid)%nvid(      ishft + jv) = zval2d(jv)
               iom_file(kiomid)%ndims(     ishft + jv) = 2
               iom_file(kiomid)%dimsz(1:2, ishft + jv) = (/ inx, iny /)
               iom_file(kiomid)%ofs(       ishft + jv) = 0.
            END DO
            ! 3d variable
            ishft = in0d + in1d + in2d 
            DO jv = 1, in3d
               iom_file(kiomid)%cn_var(    ishft + jv) = TRIM(clna3d(jv))
               iom_file(kiomid)%nvid(      ishft + jv) = zval3d(jv)
               iom_file(kiomid)%ndims(     ishft + jv) = 3
               iom_file(kiomid)%dimsz(1:3, ishft + jv) = (/ inx, iny, jpk /)
               iom_file(kiomid)%ofs(       ishft + jv) = 0.
            END DO
         ELSE                 ! new file
            iom_file(kiomid)%nvars = 0
            iom_file(kiomid)%irec  = 2
            ! store file informations
            WRITE( idrst, REC = 1, IOSTAT = ios, ERR = 987 ) irecl8, kdompar(:,1), jpk   ! store domain size
         ENDIF
      ENDIF
987   CONTINUE
      IF( ios /= 0 ) THEN
         WRITE(ctmp1,*) '           iostat = ', ios
         CALL ctl_stop( TRIM(clinfo), '   error in opening file '//TRIM(cdname), ctmp1 )
      ENDIF
      !
   END SUBROUTINE iom_rstdimg_open


   SUBROUTINE iom_rstdimg_close( kiomid )
      !!--------------------------------------------------------------------
      !!                   ***  SUBROUTINE  iom_rstdimg_close  ***
      !!
      !! ** Purpose : close an input file
      !!--------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kiomid   ! iom identifier of the file to be closed
      !
      CHARACTER(LEN=100)                      ::   clinfo                     ! info character
      INTEGER                                 ::   jv                         ! loop counter
      INTEGER                                 ::   irecl8                     ! record length
      INTEGER                                 ::   ios                        ! IO status
      INTEGER                                 ::   irhd                       ! record of the header infos
      INTEGER                                 ::   ivnum                      ! number of variables
      INTEGER                                 ::   idrst                      ! file logical unit
      INTEGER                                 ::   inx, iny, inz              ! x,y,z dimension of the variable
      INTEGER                                 ::   in0d, in1d, in2d, in3d     ! number of 0/1/2/3D variables
      CHARACTER(LEN=jpvnl), DIMENSION(jpmax_vars) ::   clna0d, clna1d, clna2d, clna3d    ! name of 0/1/2/3D variables
      REAL(wp),          DIMENSION(jpmax_vars) ::   zval0d, zval1d, zval2d, zval3d    ! value of 0d variables or record
      !                                                                               ! position for 1/2/3D variables
      !---------------------------------------------------------------------
      !
      clinfo = '                    iom_rstdimg_close ~~~  '
      idrst = iom_file(kiomid)%nfid   ! get back the logical unit of the restart file
      ! test if we can write in the file (test with INQUIRE gives alsways YES even with read only files...)
      READ(  idrst, REC = 1, IOSTAT = ios, ERR = 987 ) irecl8, inx, iny, inz   
      WRITE( idrst, REC = 1, IOSTAT = ios            ) irecl8, inx, iny, inz   
      ! We can write in the file => we update its header before closing
      IF( ios == 0 ) THEN
         READ( idrst, REC = 1, IOSTAT = ios, ERR = 987 ) irecl8, inx, iny, inz   ! get back domain size
         irhd = iom_file(kiomid)%irec
         ivnum = iom_file(kiomid)%nvars
         in0d = 0   ;   in1d = 0   ;   in2d = 0   ;   in3d = 0
         DO jv = 1, ivnum   ! loop on each variable to get its name and value/record position
            SELECT CASE (iom_file(kiomid)%ndims(jv))
            CASE (0)   ! scalar variable
               in0d = in0d + 1
               clna0d(in0d) = TRIM(iom_file(kiomid)%cn_var(jv))
               zval0d(in0d) = iom_file(kiomid)%ofs(jv)   ! warning: trick... we use ofs to store the value
            CASE (1)   ! 1d variable
               in1d = in1d + 1
               clna1d(in1d) = TRIM(iom_file(kiomid)%cn_var(jv))
               zval1d(in1d) = iom_file(kiomid)%nvid(jv)
            CASE (2)   ! 2d variable
               in2d = in2d + 1
               clna2d(in2d) = TRIM(iom_file(kiomid)%cn_var(jv))
               zval2d(in2d) = iom_file(kiomid)%nvid(jv)
            CASE (3)   ! 3d variable
               in3d = in3d + 1
               clna3d(in3d) = TRIM(iom_file(kiomid)%cn_var(jv))
               zval3d(in3d) = iom_file(kiomid)%nvid(jv)
            CASE DEFAULT   ;   CALL ctl_stop( TRIM(clinfo), 'Should not ne there...' )
            END SELECT
         END DO
         ! force to have at least 1 variable in each list (not necessary (?), but safer...)
         IF( in0d == 0 ) THEN   ;   in0d = 1   ;   clna0d(1) = 'no0d'   ;   zval0d(1) = -1.   ;   ENDIF
         IF( in1d == 0 ) THEN   ;   in1d = 1   ;   clna1d(1) = 'no1d'   ;   zval1d(1) = -1.   ;   ENDIF
         IF( in2d == 0 ) THEN   ;   in2d = 1   ;   clna2d(1) = 'no2d'   ;   zval2d(1) = -1.   ;   ENDIF
         IF( in3d == 0 ) THEN   ;   in3d = 1   ;   clna3d(1) = 'no3d'   ;   zval3d(1) = -1.   ;   ENDIF
         ! update the file header before closing it
         WRITE( idrst, REC = 1, IOSTAT = ios, ERR = 987 )              &
            &   irecl8, inx, iny, inz, in0d, in1d, in2d, in3d, -irhd,  &   ! to indicate change in dimg format, use negative irhd
            &   jpni, jpnj, jpnij, narea, jpiglo, jpjglo
         IF( (ivnum * (jpvnl + wp)) > irecl8 ) THEN 
            CALL ctl_stop( TRIM(clinfo),   &
                 &   'Last record size is too big... You could reduce the value of jpvnl' )
         ELSE 
            WRITE( idrst, REC = irhd, IOSTAT = ios, ERR = 987 )                        &
                 &   clna0d(1:in0d), zval0d(1:in0d), clna1d(1:in1d), zval1d(1:in1d),   &
                 &   clna2d(1:in2d), zval2d(1:in2d), clna3d(1:in3d), zval3d(1:in3d)
            ! save now layout array at the end of the file, 1 per record
            WRITE ( idrst, REC = irhd + 1 , IOSTAT = ios, ERR = 987 )  nlcit
            WRITE ( idrst, REC = irhd + 2 , IOSTAT = ios, ERR = 987 )  nlcjt
            WRITE ( idrst, REC = irhd + 3 , IOSTAT = ios, ERR = 987 )  nldit
            WRITE ( idrst, REC = irhd + 4 , IOSTAT = ios, ERR = 987 )  nldjt
            WRITE ( idrst, REC = irhd + 5 , IOSTAT = ios, ERR = 987 )  nleit
            WRITE ( idrst, REC = irhd + 6 , IOSTAT = ios, ERR = 987 )  nlejt
            WRITE ( idrst, REC = irhd + 7 , IOSTAT = ios, ERR = 987 )  nimppt
            WRITE ( idrst, REC = irhd + 8 , IOSTAT = ios, ERR = 987 )  njmppt
         ENDIF
      ELSE
         ios = 0   ! we cannot write in the file
      ENDIF
      !
      CLOSE( idrst, IOSTAT = ios, ERR = 987 )
987   CONTINUE
      IF( ios /= 0 ) THEN
         WRITE(ctmp1,*) '           iostat = ', ios
         CALL ctl_stop( TRIM(clinfo),   &
            &   '   error when updating the header of '//TRIM(iom_file(kiomid)%name), ctmp1 )
      ENDIF
      !    
   END SUBROUTINE iom_rstdimg_close


   SUBROUTINE iom_rstdimg_g0d( kiomid, kvid, pvar )
      !!-----------------------------------------------------------------------
      !!                  ***  ROUTINE  iom_rstdimg_g0d  ***
      !!
      !! ** Purpose : read a scalar with RSTDIMG
      !!-----------------------------------------------------------------------
      INTEGER,  INTENT(in   ) ::   kiomid    ! Identifier of the file
      INTEGER,  INTENT(in   ) ::   kvid      ! variable id
      REAL(wp), INTENT(  out) ::   pvar      ! read field
      !---------------------------------------------------------------------
      !
      pvar = iom_file(kiomid)%ofs(kvid)   ! warning: trick... we use ofs to store the value
      !
   END SUBROUTINE iom_rstdimg_g0d


   SUBROUTINE iom_rstdimg_rp0d( kiomid, cdvar, kvid, pv_r0d )
      !!--------------------------------------------------------------------
      !!                   ***  SUBROUTINE  iom_rstdimg_rstput  ***
      !!
      !! ** Purpose : write a scalar with RSTDIMG
      !!--------------------------------------------------------------------
      INTEGER                   , INTENT(in) ::   kiomid   ! Identifier of the file 
      CHARACTER(len=*)          , INTENT(in) ::   cdvar    ! time axis name
      INTEGER                   , INTENT(in) ::   kvid     ! variable id
      REAL(wp)                  , INTENT(in) ::   pv_r0d   ! written 0d field
      !
      CHARACTER(LEN=100) ::   clinfo     ! info character
      INTEGER            ::   idvar      ! variable id
      !---------------------------------------------------------------------
      !  
      clinfo = '                    iom_rstdimg_rp0d ~~~  '
      IF( kvid <= 0 ) THEN   !   new variable
         idvar = iom_file(kiomid)%nvars + 1
      ELSE                   !   the variable already exists in the file
         idvar = kvid
      ENDIF
      IF( idvar <= jpmax_vars ) THEN
         iom_file(kiomid)%nvars = idvar
         iom_file(kiomid)%cn_var(idvar) = TRIM(cdvar)
         iom_file(kiomid)%nvid(   idvar) = 1   ! useless, Od variables a strored in record 1
         iom_file(kiomid)%ndims(  idvar) = 0
         iom_file(kiomid)%luld(   idvar) = .FALSE.
         iom_file(kiomid)%scf(    idvar) = 1.
         iom_file(kiomid)%ofs(    idvar) = pv_r0d   ! warning: trick... we use ofs to store the value
      ELSE
         CALL ctl_stop( TRIM(clinfo), 'increase the value of jpmax_vars' )
      ENDIF
   END SUBROUTINE iom_rstdimg_rp0d


   SUBROUTINE iom_rstdimg_g123d( kiomid, kdom  , kvid, kx1, kx2, ky1, ky2,   &
         &                       pv_r1d, pv_r2d, pv_r3d )
      !!-----------------------------------------------------------------------
      !!                  ***  ROUTINE  iom_rstdimg_g123d  ***
      !!
      !! ** Purpose : read a 1D/2D/3D variable with RSTDIMG
      !!
      !! ** Method : read ONE record at each CALL
      !!-----------------------------------------------------------------------
      INTEGER                    , INTENT(in   )           ::   kiomid     ! iom identifier of the file
      INTEGER                    , INTENT(in   )           ::   kdom       ! Type of domain to be read
      INTEGER                    , INTENT(in   )           ::   kvid       ! variable id
      INTEGER ,                    INTENT(inout)           ::   kx1, kx2, ky1, ky2   ! subdomain indexes
      REAL(wp), DIMENSION(:)     , INTENT(  out), OPTIONAL ::   pv_r1d     ! read field (1D case)
      REAL(wp), DIMENSION(:,:)   , INTENT(  out), OPTIONAL ::   pv_r2d     ! read field (2D case)
      REAL(wp), DIMENSION(:,:,:) , INTENT(  out), OPTIONAL ::   pv_r3d     ! read field (3D case)

      CHARACTER(LEN=100) ::   clinfo               ! info character
      INTEGER            ::   ios                  ! IO status
      INTEGER            ::   jk                   ! loop counter
      INTEGER            ::   idrst                ! logical unit of the restart file
      !---------------------------------------------------------------------
      clinfo = '                    iom_rstdimg_g123d ~~~  '
      !
!JMM      IF( kdom == jpdom_data .OR. kdom == jpdom_global ) THEN
!JMM         CALL ctl_stop( TRIM(clinfo), TRIM(iom_file(kiomid)%cn_var(kvid))//': case not coded for rstdimg files' )
!JMM      ELSE
      !
         idrst = iom_file(kiomid)%nfid   ! get back the logical unit of the restart file
         ! modify the subdomain indexes because we cannot directly extract the appropriate subdomaine
         IF(     kdom == jpdom_local_full    ) THEN   ;   kx1 = 1   ;   kx2 = jpi    ;   ky1 = 1
         ELSEIF( kdom == jpdom_local_noextra ) THEN   ;   kx1 = 1   ;   kx2 = nlci   ;   ky1 = 1
!JMM {
         ELSEIF( kdom == jpdom_data          ) THEN   ;   kx1 = 1   ;   kx2 = jpidta ;   ky1 = 1 ; ky2 = jpjdta
         ELSEIF( kdom == jpdom_global        ) THEN   ;   kx1 = 1   ;   kx2 = jpiglo ;   ky1 = 1 ; ky2 = jpjglo
!JMM }
         ENDIF
         !
         IF(     PRESENT(pv_r1d) ) THEN   ! read 1D variables
            READ(    idrst, REC = iom_file(kiomid)%nvid(kvid)         , IOSTAT = ios, ERR = 987 )   pv_r1d(:)
         ELSEIF( PRESENT(pv_r2d) ) THEN   ! read 2D variables
            READ(    idrst, REC = iom_file(kiomid)%nvid(kvid)         , IOSTAT = ios, ERR = 987 )   pv_r2d(kx1:kx2, ky1:ky2    )
         ELSEIF( PRESENT(pv_r3d) ) THEN   ! read 3D variables
            DO jk = 1, iom_file(kiomid)%dimsz(3,kvid)   ! do loop on each level
               READ( idrst, REC = iom_file(kiomid)%nvid(kvid) + jk - 1, IOSTAT = ios, ERR = 987 )   pv_r3d(kx1:kx2, ky1:ky2, jk)
            END DO
         ENDIF
987      CONTINUE
         IF( ios /= 0 ) THEN
            WRITE(ctmp1,*) '           iostat = ', ios
            CALL ctl_stop( TRIM(clinfo), '   IO error with file '//TRIM(iom_file(kiomid)%name), ctmp1 )
         ENDIF
!JMM      ENDIF
      !
   END SUBROUTINE iom_rstdimg_g123d


   SUBROUTINE iom_rstdimg_rp123d( kiomid, cdvar, kvid, pv_r1d, pv_r2d, pv_r3d )
      !!--------------------------------------------------------------------
      !!                   ***  SUBROUTINE  iom_rstdimg_rstput  ***
      !!
      !! ** Purpose : write a 2D/3D variable with RSTDIMG
      !!--------------------------------------------------------------------
      INTEGER                         , INTENT(in)           ::   kiomid   ! Identifier of the file 
      CHARACTER(len=*)                , INTENT(in)           ::   cdvar    ! time axis name
      INTEGER                         , INTENT(in)           ::   kvid     ! variable id
      REAL(wp), DIMENSION(          :), INTENT(in), OPTIONAL ::   pv_r1d   ! written 1d field
      REAL(wp), DIMENSION(:  ,:      ), INTENT(in), OPTIONAL ::   pv_r2d   ! written 2d field
      REAL(wp), DIMENSION(:  ,:  ,:  ), INTENT(in), OPTIONAL ::   pv_r3d   ! written 3d field
      !
      CHARACTER(LEN=100)    ::   clinfo               ! info character
      INTEGER               ::   irecl8               ! reacord length
      INTEGER               ::   ios                  ! IO status
      INTEGER               ::   idrst                ! reacord length
      INTEGER               ::   inx, iny, inz        ! x,y,z dimension of the variable
      INTEGER               ::   idvar                ! variable id
      INTEGER               ::   istop                ! temporary storage of nstop
      INTEGER               ::   irec                 ! record number
      INTEGER               ::   ix1, ix2, iy1, iy2   ! subdomain indexes
      INTEGER               ::   jk                   ! loop counter
      !---------------------------------------------------------------------
      !
      clinfo = '          iom_rstdimg_rp123d, file: '//TRIM(iom_file(kiomid)%name)//', var: '//TRIM(cdvar)
      istop = nstop                   ! store the actual value of nstop
      irec = iom_file(kiomid)%irec    ! get back the record number of the variable
      idrst = iom_file(kiomid)%nfid   ! get back the logical unit of the restart file
      IF( kvid <= 0 ) THEN   !   new variable
         idvar = iom_file(kiomid)%nvars + 1
      ELSE                   !   the variable already exists in the file
         idvar = kvid
      ENDIF
      IF( idvar > jpmax_vars )   CALL ctl_stop( TRIM(clinfo), 'increase the value of jpmax_vars' )
      IF( .NOT. PRESENT(pv_r1d) ) THEN
         ! find which part of data must be written
         READ( idrst, REC = 1, IOSTAT = ios, ERR = 987 ) irecl8, inx, iny, inz
         IF(     inx == (nlei - nldi + 1) .AND. iny == (nlej - nldj + 1) ) THEN
            ix1 = nldi   ;   ix2 = nlei   ;   iy1 = nldj   ;   iy2 = nlej
         ELSEIF( inx == nlci              .AND. iny == nlcj              ) THEN
            ix1 = 1      ;   ix2 = nlci   ;   iy1 = 1      ;   iy2 = nlcj
         ELSEIF( inx == jpi               .AND. iny == jpj               ) THEN
            ix1 = 1      ;   ix2 = jpi    ;   iy1 = 1      ;   iy2 = jpj
         ELSE 
            CALL ctl_stop( clinfo, 'should have been an impossible case...' )
         ENDIF
      ENDIF
      IF( istop == nstop ) THEN
         ! write the data
         IF(     PRESENT(pv_r1d) ) THEN   ! 1D variable
            WRITE( idrst, REC = irec            , IOSTAT = ios, ERR = 987 ) pv_r1d(:)
         ELSEIF( PRESENT(pv_r2d) ) THEN   ! 2D variable
            WRITE( idrst, REC = irec            , IOSTAT = ios, ERR = 987 ) pv_r2d(ix1:ix2, iy1:iy2    )
         ELSEIF( PRESENT(pv_r3d) ) THEN   ! 3D variable
            DO jk = 1, jpk   ! do loop on each level
               WRITE( idrst, REC = irec + jk - 1, IOSTAT = ios, ERR = 987 ) pv_r3d(ix1:ix2, iy1:iy2, jk)
            END DO
         ENDIF
         ! fill the file informations
         iom_file(kiomid)%nvars = idvar
         IF(     PRESENT(pv_r1d) ) THEN
            iom_file(kiomid)%irec              = irec + 1
            iom_file(kiomid)%ndims(     idvar) = 1
            iom_file(kiomid)%dimsz(1  , idvar) = inz
         ELSEIF( PRESENT(pv_r2d) ) THEN
            iom_file(kiomid)%irec              = irec + 1
            iom_file(kiomid)%ndims(     idvar) = 2
            iom_file(kiomid)%dimsz(1:2, idvar) = (/ inx, iny /)
         ELSEIF( PRESENT(pv_r3d) ) THEN
            iom_file(kiomid)%irec              = irec + inz
            iom_file(kiomid)%ndims(     idvar) = 3
            iom_file(kiomid)%dimsz(1:3, idvar) = (/ inx, iny, inz /)
         ENDIF
         iom_file(kiomid)%cn_var(idvar) = TRIM(cdvar)
         iom_file(kiomid)%nvid(  idvar) = irec
         iom_file(kiomid)%luld(  idvar) = .FALSE.
         iom_file(kiomid)%scf(   idvar) = 1.
         iom_file(kiomid)%ofs(   idvar) = 0.
      ENDIF
987   CONTINUE
      IF( ios /= 0 ) THEN
         WRITE(ctmp1,*) '           iostat = ', ios
         CALL ctl_stop( TRIM(clinfo), '   IO error with file '//TRIM(iom_file(kiomid)%name), ctmp1 )
      ELSE
         IF(lwp) WRITE(numout,*) TRIM(clinfo)//' written ok'
      ENDIF
      !     
   END SUBROUTINE iom_rstdimg_rp123d


   !!======================================================================
END MODULE iom_rstdimg
