MODULE diadimg
   !!======================================================================
   !!                     ***  MODULE  diadimg  ***
   !! Ocean diagnostics :  write ocean output files in dimg direct access format (mpp)
   !!=====================================================================
# if defined key_dimgout
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers 
   USE dom_oce         ! ocean space and time domain
   USE in_out_manager  ! I/O manager
   USE daymod          ! calendar
   USE lib_mpp

   IMPLICIT NONE
   PRIVATE

   PUBLIC dia_wri_dimg            ! called by trd_mld (eg)
   PUBLIC dia_wri_dimg_alloc      ! called by nemo_alloc in nemogcm.F90


   !! These workspace arrays are inside the module so that we can make them
   !! allocatable in a clean way. Not done in wrk_nemo because these are of KIND(sp).
   REAL(sp), ALLOCATABLE, SAVE, DIMENSION(:,:) :: z42d    ! 2d temporary workspace (sp)
   REAL(sp), ALLOCATABLE, SAVE, DIMENSION(:)   :: z4dep   ! vertical level (sp)

   !! * Substitutions
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: diadimg.F90 4292 2013-11-20 16:28:04Z cetlod $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   FUNCTION dia_wri_dimg_alloc()
      !!---------------------------------------------------------------------
      !!        *** ROUTINE dia_wri_dimg_alloc ***
      !!
      !!---------------------------------------------------------------------
      INTEGER :: dia_wri_dimg_alloc   ! return value
      !!---------------------------------------------------------------------
      !
      IF( .NOT. ALLOCATED( z42d ) )THEN

         ALLOCATE( z42d(jpi,jpj), z4dep(jpk), STAT=dia_wri_dimg_alloc )

         IF( lk_mpp                  )   CALL mpp_sum ( dia_wri_dimg_alloc )
         IF( dia_wri_dimg_alloc /= 0 )   CALL ctl_warn('dia_wri_dimg_alloc: allocation of array failed.')

      ELSE

         dia_wri_dimg_alloc = 0

      ENDIF
      !
  END FUNCTION dia_wri_dimg_alloc


  SUBROUTINE dia_wri_dimg(cd_name, cd_text, ptab, klev, cd_type , ksubi )
    !!-------------------------------------------------------------------------
    !!        *** ROUTINE dia_wri_dimg ***
    !!
    !! ** Purpose : write ptab in the dimg file cd_name, with comment cd_text.
    !!       ptab has klev x 2D fields
    !!
    !! ** Action :   Define header variables from the config parameters
    !!       Open the dimg file on unit inum = 14 ( IEEE I4R4 file )
    !!       Write header on record 1
    !!       Write ptab on the following klev records
    !!
    !! History :  2003-12 (J.M. Molines ) : Original. Replace ctl_opn, writn2d
    !!---------------------------------------------------------------------------
    CHARACTER(len=*),INTENT(in) ::   &
         &                            cd_name,  &  ! dimg file name
         &                            cd_text      ! comment to write on record #1
    INTEGER, INTENT(in) ::            klev         ! number of level in ptab to write
    REAL(wp),INTENT(in), DIMENSION(:,:,:) :: ptab  ! 3D array to write 
    CHARACTER(LEN=1),INTENT(in) ::    cd_type      ! either 'T', 'W' or '2' , depending on the vertical
    !                                              ! grid for ptab. 2 stands for 2D file
    INTEGER, INTENT(in), OPTIONAL, DIMENSION(klev) :: ksubi 

    !! * Local declarations
    INTEGER :: jk, jn           ! dummy loop indices
    INTEGER :: irecl4,             &    ! record length in bytes
         &       inum,             &    ! logical unit (set to 14)
         &       irec,             &    ! current record to be written
         &       irecend                ! record number where nclit... are stored
    REAL(sp)                    :: zdx,zdy,zspval,zwest,ztimm
    REAL(sp)                    :: zsouth, zdtwr

    CHARACTER(LEN=80) :: clname                ! name of file in case of dimgnnn
    CHARACTER(LEN=4) :: clver='@!01'           ! dimg string identifier
    !!---------------------------------------------------------------------------

    !                                      ! allocate dia_wri_dimg array
    IF( dia_wri_dimg_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'dia_wri_dimg : unable to allocate arrays' )

    !! * Initialisations

    !{ DRAKKAR : add jpnij*sp in the test because for many procs, this can be the limit
    irecl4 = MAX(jpi*jpj*sp , 84+(18+1+jpk)*sp, jpnij*sp )
!!!! jm34?   inum = 14
    !}

    zspval=0.0_sp    ! special values on land
    !  the 'numerical' grid is described. The geographical one is in a grid file
    zdx=1._sp
    zdy=1._sp
    zsouth=njmpp * 1._sp
    zwest=nimpp * 1._sp
    zdtwr=nn_write * rn_rdt  / 86400_wp   ! write step in days
    !  time in days since the historical begining of the run (nit000 = 0 ) 
    ztimm=adatrj - zdtwr/2.
    

    SELECT CASE ( cd_type )

    CASE ( 'T')
       z4dep(:)=gdept_1d(:)

    CASE ( 'W' )
       z4dep(:)=gdepw_1d(:)

    CASE ( '2' )
       z4dep(1:klev) =(/(jk, jk=1,klev)/)

    CASE ( 'I' )
       z4dep(1:klev) = ksubi(1:klev)

    CASE DEFAULT
       IF(lwp) WRITE(numout,*) ' E R R O R : bad cd_type in dia_wri_dimg '
       STOP 'dia_wri_dimg'

    END SELECT

    IF ( ln_dimgnnn  ) THEN
     !{ DRAKKAR : fix when writing Daily files
     clver='@!01'           ! dimg string identifier
     !}
     irecl4 = MAX(jpi*jpj*sp , 84+(18+jpk)*sp + 8*jpnij*sp +3 *sp  )
       WRITE(clname,'(a,a,i3.3)') TRIM(cd_name),'.',narea
       CALL ctl_opn(inum, clname,'UNKNOWN','UNFORMATTED','DIRECT',irecl4,numout,lwp, cdirout=cn_dirout)
       WRITE(inum,REC=1 ) clver, cd_text, irecl4, &
            &     jpi,jpj, klev, 1 , 1 ,            &
            &     zwest, zsouth, zdx, zdy, zspval,  &
            &     z4dep(1:klev),                    &
            &     ztimm,                            &
            &     narea, jpnij,jpiglo,jpjglo,jpizoom, jpjzoom,    &    ! extension to dimg for mpp output
            &     nlcit,nlcjt, nldit, nldjt, nleit, nlejt, nimppt, njmppt, &   !
            &     ndate0, zdtwr, nleapy

       !! * Write klev levels
       IF ( cd_type == 'I' ) THEN

          DO jk = 1, klev
             irec =1 + jk
             z42d(:,:) = ptab(:,:,ksubi(jk))
             WRITE(inum,REC=irec)  z42d(:,:)
          END DO
       ELSE
          DO jk = 1, klev
             irec =1 + jk
             z42d(:,:) = ptab(:,:,jk)
             WRITE(inum,REC=irec)  z42d(:,:)
          END DO
       ENDIF
    ELSE
       clver='@!04'           ! dimg string identifier
       ! note that version @!02 is optimized with respect to record length.
       ! The vertical dep variable is reduced to klev instead of klev*jpnij :
       !   this is OK for jpnij < 181 (jpk=46)
       ! for more processors, irecl4 get huge and that's why we switch to '@!03':
       !  In this case we just add an extra integer to the standard dimg structure,
       !  which is a record number where the arrays nlci etc... starts (1 per record)
       ! version '@!04 add an extra record with ndate0 and nn_write*rn_rdt
       
       !! Standard dimgproc (1 file per variable, all procs. write to this file )
       !! * Open file
       CALL ctl_opn(inum, cd_name,'UNKNOWN','UNFORMATTED','DIRECT',irecl4,numout,lwp, cdirout=cn_dirout)

       !! * Write header on record #1
       irecend=1 + klev*jpnij 
       IF(lwp) WRITE(inum,REC=1 ) clver, cd_text, irecl4, &
            &     jpi,jpj, klev, 1 , 1 ,            &
            &     zwest, zsouth, zdx, zdy, zspval,  &
            &     z4dep(1:klev),       &
            &     ztimm,                            &
            &     narea, jpnij,jpiglo,jpjglo,jpizoom, jpjzoom, irecend 
       IF (lwp ) THEN
         WRITE(inum,REC=irecend + 1 ) nlcit
         WRITE(inum,REC=irecend + 2 ) nlcjt
         WRITE(inum,REC=irecend + 3 ) nldit
         WRITE(inum,REC=irecend + 4 ) nldjt
         WRITE(inum,REC=irecend + 5 ) nleit
         WRITE(inum,REC=irecend + 6 ) nlejt
         WRITE(inum,REC=irecend + 7 ) nimppt
         WRITE(inum,REC=irecend + 8 ) njmppt
         WRITE(inum,REC=irecend + 9 ) ndate0, zdtwr, nleapy
       ENDIF
      !   &    ! extension to dimg for mpp output
      !   &     nlcit,nlcjt, nldit, nldjt, nleit, nlejt, nimppt, njmppt  !

       !! * Write klev levels
       IF ( cd_type == 'I' ) THEN

          DO jk = 1, klev
             irec =1 + klev * (narea -1) + jk
             z42d(:,:) = ptab(:,:,ksubi(jk))
             WRITE(inum,REC=irec)  z42d(:,:)
          END DO
       ELSE
          DO jk = 1, klev
             irec =1 + klev * (narea -1) + jk
             z42d(:,:) = ptab(:,:,jk)
             WRITE(inum,REC=irec)  z42d(:,:)
          END DO
       ENDIF
    ENDIF

    !! * Close the file
    CLOSE(inum)

  END SUBROUTINE dia_wri_dimg

#  else
   !!----------------------------------------------------------------------
   !!   Default option :                                       Empty module
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dia_wri_dimg( cd_name, cd_exper, ptab, klev, cd_type )
      REAL, DIMENSION(:,:,:) :: ptab
      INTEGER :: klev
      CHARACTER(LEN=80) :: cd_name, cd_exper,cd_type
      WRITE(*,*) ' This print must never occur ', cd_name, cd_exper,ptab, klev, cd_type
      WRITE(*,*) ' this routine is here just for compilation '
   END SUBROUTINE dia_wri_dimg
# endif
   !!======================================================================
END MODULE diadimg
