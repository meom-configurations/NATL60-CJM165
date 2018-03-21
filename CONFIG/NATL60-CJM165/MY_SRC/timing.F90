MODULE timing
   !!========================================================================
   !!                     ***  MODULE  timing  ***
   !!========================================================================
   !! History : 4.0  ! 2001-05  (R. Benshila)   
   !!------------------------------------------------------------------------

   !!------------------------------------------------------------------------
   !!   timming_init    : initialize timing process 
   !!   timing_start    : start Timer
   !!   timing_stop     : stop  Timer
   !!   timing_reset    : end timing variable creation
   !!   timing_finalize : compute stats and write output in calling w*_info 
   !!   timing_ini_var  : create timing variables 
   !!   timing_listing  : print instumented subroutines in ocean.output
   !!   wcurrent_info   : compute and print detailed stats on the current CPU
   !!   wave_info       : compute and print averaged statson all processors
   !!   wmpi_info       : compute and write global stats  
   !!   supress         : suppress an element of the timing linked list  
   !!   insert          : insert an element of the timing linked list  
   !!------------------------------------------------------------------------
   USE in_out_manager  ! I/O manager 
   USE dom_oce         ! ocean domain
   USE lib_mpp          
   
   IMPLICIT NONE
   PRIVATE

   PUBLIC   timing_init, timing_finalize   ! called in nemogcm module 
   PUBLIC   timing_reset                   ! called in step module 
   PUBLIC   timing_start, timing_stop      ! called in each routine to time 
   
#if defined key_mpp_mpi
   INCLUDE 'mpif.h'
#endif

   ! Variables for fine grain timing
   TYPE timer
      CHARACTER(LEN=20)  :: cname
   	  REAL(wp)  :: t_cpu, t_clock, tsum_cpu, tsum_clock, tmax_cpu, tmax_clock, tmin_cpu, tmin_clock, tsub_cpu, tsub_clock
      INTEGER :: ncount, ncount_max, ncount_rate  
      INTEGER :: niter
      LOGICAL :: l_tdone
      TYPE(timer), POINTER :: next => NULL()
      TYPE(timer), POINTER :: prev => NULL()
      TYPE(timer), POINTER :: parent_section => NULL()
   END TYPE timer
    
   TYPE alltimer
      CHARACTER(LEN=20), DIMENSION(:), POINTER :: cname => NULL()
   	  REAL(wp), DIMENSION(:), POINTER :: tsum_cpu   => NULL()
   	  REAL(wp), DIMENSION(:), POINTER :: tsum_clock => NULL()
   	  INTEGER, DIMENSION(:), POINTER :: niter => NULL()
      TYPE(alltimer), POINTER :: next => NULL()
      TYPE(alltimer), POINTER :: prev => NULL()
   END TYPE alltimer 
 
   TYPE(timer), POINTER :: s_timer_root => NULL()
   TYPE(timer), POINTER :: s_timer      => NULL()
   TYPE(timer), POINTER :: s_wrk        => NULL()
   REAL(wp) :: t_overclock, t_overcpu
   LOGICAL :: l_initdone = .FALSE.
   INTEGER :: nsize
   
   ! Variables for coarse grain timing
   REAL(wp) :: tot_etime, tot_ctime
   REAL(kind=wp), DIMENSION(2)     :: t_elaps, t_cpu
   REAL(wp), ALLOCATABLE, DIMENSION(:) :: all_etime, all_ctime
   INTEGER :: nfinal_count, ncount, ncount_rate, ncount_max
   INTEGER, DIMENSION(8)           :: nvalues
   CHARACTER(LEN=8), DIMENSION(2)  :: cdate
   CHARACTER(LEN=10), DIMENSION(2) :: ctime
   CHARACTER(LEN=5)                :: czone
    
   ! From of ouput file (1/proc or one global)   !RB to put in nammpp or namctl
   LOGICAL :: ln_onefile = .TRUE. 
   LOGICAL :: lwriter
   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO Consortium (2011)
   !! $Id: timing.F90 5120 2015-03-03 16:11:55Z acc $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE timing_start(cdinfo)
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE timing_start  ***
      !! ** Purpose :   collect execution time
      !!----------------------------------------------------------------------
      CHARACTER(len=*), INTENT(in) :: cdinfo
      !
       
      ! Create timing structure at first call
      IF( .NOT. l_initdone ) THEN
         CALL timing_ini_var(cdinfo)
      ELSE
         s_timer => s_timer_root
         DO WHILE( TRIM(s_timer%cname) /= TRIM(cdinfo) ) 
            IF( ASSOCIATED(s_timer%next) ) s_timer => s_timer%next
         END DO
      ENDIF         
      s_timer%l_tdone = .FALSE.
      s_timer%niter = s_timer%niter + 1
      s_timer%t_cpu = 0.
      s_timer%t_clock = 0.
                  
      ! CPU time collection
      CALL CPU_TIME( s_timer%t_cpu  )
      ! clock time collection
#if defined key_mpp_mpi
      s_timer%t_clock= MPI_Wtime()
#else
      CALL SYSTEM_CLOCK(COUNT_RATE=s_timer%ncount_rate, COUNT_MAX=s_timer%ncount_max)
      CALL SYSTEM_CLOCK(COUNT = s_timer%ncount)
#endif
      !
   END SUBROUTINE timing_start


   SUBROUTINE timing_stop(cdinfo, csection)
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE timing_stop  ***
      !! ** Purpose :   finalize timing and output
      !!----------------------------------------------------------------------
      CHARACTER(len=*), INTENT(in) :: cdinfo
      CHARACTER(len=*), INTENT(in), OPTIONAL :: csection
      !
      INTEGER  :: ifinal_count, iperiods    
      REAL(wp) :: zcpu_end, zmpitime
      !
      s_wrk => NULL()

      ! clock time collection
#if defined key_mpp_mpi
      zmpitime = MPI_Wtime()
#else
      CALL SYSTEM_CLOCK(COUNT = ifinal_count)
#endif
      ! CPU time collection
      CALL CPU_TIME( zcpu_end )

      s_timer => s_timer_root
      DO WHILE( TRIM(s_timer%cname) /= TRIM(cdinfo) ) 
         IF( ASSOCIATED(s_timer%next) ) s_timer => s_timer%next
      END DO
 
      ! CPU time correction
      s_timer%t_cpu  = zcpu_end - s_timer%t_cpu - t_overcpu - s_timer%tsub_cpu
  
      ! clock time correction
#if defined key_mpp_mpi
      s_timer%t_clock = zmpitime - s_timer%t_clock - t_overclock - s_timer%tsub_clock
#else
      iperiods = ifinal_count - s_timer%ncount
      IF( ifinal_count < s_timer%ncount )  &
          iperiods = iperiods + s_timer%ncount_max 
      s_timer%t_clock  = REAL(iperiods) / s_timer%ncount_rate - t_overclock - s_timer%tsub_clock
#endif
      
      ! Correction of parent section
      IF( .NOT. PRESENT(csection) ) THEN
         s_wrk => s_timer
         DO WHILE ( ASSOCIATED(s_wrk%parent_section ) )
            s_wrk => s_wrk%parent_section
            s_wrk%tsub_cpu   = s_wrk%tsub_cpu   + s_timer%t_cpu 
            s_wrk%tsub_clock = s_wrk%tsub_clock + s_timer%t_clock              
         END DO
      ENDIF
            
      ! time diagnostics 
      s_timer%tsum_clock = s_timer%tsum_clock + s_timer%t_clock 
      s_timer%tsum_cpu   = s_timer%tsum_cpu   + s_timer%t_cpu
!RB to use to get min/max during a time integration
!      IF( .NOT. l_initdone ) THEN
!         s_timer%tmin_clock = s_timer%t_clock 
!         s_timer%tmin_cpu   = s_timer%t_cpu 
!      ELSE
!         s_timer%tmin_clock = MIN( s_timer%tmin_clock, s_timer%t_clock ) 
!         s_timer%tmin_cpu   = MIN( s_timer%tmin_cpu  , s_timer%t_cpu   ) 
!      ENDIF   
!      s_timer%tmax_clock = MAX( s_timer%tmax_clock, s_timer%t_clock ) 
!      s_timer%tmax_cpu   = MAX( s_timer%tmax_cpu  , s_timer%t_cpu   )  
      !
      s_timer%tsub_clock = 0.
      s_timer%tsub_cpu = 0.
      s_timer%l_tdone = .TRUE.
      !
   END SUBROUTINE timing_stop
 
 
   SUBROUTINE timing_init
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE timing_init  ***
      !! ** Purpose :   open timing output file
      !!----------------------------------------------------------------------
      INTEGER :: iperiods, istart_count, ifinal_count
      REAL(wp) :: zdum
      LOGICAL :: ll_f
             
      IF( ln_onefile ) THEN
         IF( lwp) CALL ctl_opn( numtime, 'timing.output', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, numout,.TRUE., narea )
         lwriter = lwp
      ELSE
         CALL ctl_opn( numtime, 'timing.output', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, numout,.FALSE., narea )
         lwriter = .TRUE.
      ENDIF
      
      IF( lwriter) THEN      
         WRITE(numtime,*)
         WRITE(numtime,*) '      CNRS - NERC - Met OFFICE - MERCATOR-ocean - CMCC - INGV'
         WRITE(numtime,*) '                             NEMO team'
         WRITE(numtime,*) '                  Ocean General Circulation Model'
         WRITE(numtime,*) '                        version 3.6  (2015) '
         WRITE(numtime,*)
         WRITE(numtime,*) '                        Timing Informations '
         WRITE(numtime,*)
         WRITE(numtime,*)
      ENDIF   
      
      ! Compute clock function overhead
#if defined key_mpp_mpi        
      t_overclock = MPI_WTIME()
      t_overclock = MPI_WTIME() - t_overclock
#else        
      CALL SYSTEM_CLOCK(COUNT_RATE=ncount_rate, COUNT_MAX=ncount_max)
      CALL SYSTEM_CLOCK(COUNT = istart_count)
      CALL SYSTEM_CLOCK(COUNT = ifinal_count)
      iperiods = ifinal_count - istart_count
      IF( ifinal_count < istart_count )  &
          iperiods = iperiods + ncount_max 
      t_overclock = REAL(iperiods) / ncount_rate
#endif

      ! Compute cpu_time function overhead
      CALL CPU_TIME(zdum)
      CALL CPU_TIME(t_overcpu)
      
      ! End overhead omputation  
      t_overcpu = t_overcpu - zdum        
      t_overclock = t_overcpu + t_overclock        

      ! Timing on date and time
      CALL DATE_AND_TIME(cdate(1),ctime(1),czone,nvalues)
    
      CALL CPU_TIME(t_cpu(1))      
#if defined key_mpp_mpi        
      ! Start elapsed and CPU time counters
      t_elaps(1) = MPI_WTIME()
#else
      CALL SYSTEM_CLOCK(COUNT_RATE=ncount_rate, COUNT_MAX=ncount_max)
      CALL SYSTEM_CLOCK(COUNT = ncount)
#endif                 
      !
   END SUBROUTINE timing_init


   SUBROUTINE timing_finalize
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE timing_finalize ***
      !! ** Purpose :  compute average time 
      !!               write timing output file
      !!----------------------------------------------------------------------
      TYPE(timer), POINTER :: s_temp
      INTEGER :: idum, iperiods, icode
      LOGICAL :: ll_ord, ll_averep
      CHARACTER(len=120) :: clfmt            
      
      ll_averep = .TRUE.
    
      ! total CPU and elapse
      CALL CPU_TIME(t_cpu(2))
      t_cpu(2)   = t_cpu(2)    - t_cpu(1)   - t_overcpu
#if defined key_mpp_mpi
      t_elaps(2) = MPI_WTIME() - t_elaps(1) - t_overclock
#else
      CALL SYSTEM_CLOCK(COUNT = nfinal_count)
      iperiods = nfinal_count - ncount
      IF( nfinal_count < ncount )  &
          iperiods = iperiods + ncount_max 
      t_elaps(2) = REAL(iperiods) / ncount_rate - t_overclock
#endif      

      ! End of timings on date & time
      CALL DATE_AND_TIME(cdate(2),ctime(2),czone,nvalues)
       
      ! Compute the numer of routines
      nsize = 0 
      s_timer => s_timer_root
      DO WHILE( ASSOCIATED(s_timer) )
         nsize = nsize + 1
         s_timer => s_timer%next
      END DO
      idum = nsize
      IF(lk_mpp) CALL mpp_sum(idum)
      IF( idum/jpnij /= nsize ) THEN
         IF( lwriter ) WRITE(numtime,*) '        ===> W A R N I N G: '
         IF( lwriter ) WRITE(numtime,*) ' Some CPU have different number of routines instrumented for timing'
         IF( lwriter ) WRITE(numtime,*) ' No detailed report on averaged timing can be provided'
         IF( lwriter ) WRITE(numtime,*) ' The following detailed report only deals with the current processor'
         IF( lwriter ) WRITE(numtime,*)
         ll_averep = .FALSE.
      ENDIF   

#if defined key_mpp_mpi      
      ! in MPI gather some info
      ALLOCATE( all_etime(jpnij), all_ctime(jpnij) )
      CALL MPI_ALLGATHER(t_elaps(2), 1, MPI_DOUBLE_PRECISION,   &
                         all_etime , 1, MPI_DOUBLE_PRECISION,   &
                         MPI_COMM_OPA, icode)
      CALL MPI_ALLGATHER(t_cpu(2) , 1, MPI_DOUBLE_PRECISION,   &
                         all_ctime, 1, MPI_DOUBLE_PRECISION,   &
                         MPI_COMM_OPA, icode)
      tot_etime = SUM(all_etime(:))
      tot_ctime = SUM(all_ctime(:))
#else
      tot_etime = t_elaps(2)
      tot_ctime = t_cpu  (2)           
#endif

      ! write output file
      IF( lwriter ) WRITE(numtime,*) 'Total timing (sum) :'
      IF( lwriter ) WRITE(numtime,*) '--------------------'
      IF( lwriter ) WRITE(numtime,"('Elapsed Time (s)  CPU Time (s)')")
      IF( lwriter ) WRITE(numtime,'(5x,f13.3,1x,f12.3)')  tot_etime, tot_ctime
      IF( lwriter ) WRITE(numtime,*) 
#if defined key_mpp_mpi
      IF( ll_averep ) CALL waver_info
      CALL wmpi_info
#endif      
      IF( lwriter ) CALL wcurrent_info
      
      clfmt='(1X,"Timing started on ",2(A2,"/"),A4," at ",2(A2,":"),A2," MET ",A3,":",A2," from GMT")'
      IF( lwriter ) WRITE(numtime, TRIM(clfmt)) &           
      &       cdate(1)(7:8), cdate(1)(5:6), cdate(1)(1:4),   &
      &       ctime(1)(1:2), ctime(1)(3:4), ctime(1)(5:6),   &
      &       czone(1:3),    czone(4:5)                     
      clfmt='(1X,  "Timing   ended on ",2(A2,"/"),A4," at ",2(A2,":"),A2," MET ",A3,":",A2," from GMT")'
      IF( lwriter ) WRITE(numtime, TRIM(clfmt)) &           
      &       cdate(2)(7:8), cdate(2)(5:6), cdate(2)(1:4),   &
      &       ctime(2)(1:2), ctime(2)(3:4), ctime(2)(5:6),   &
      &       czone(1:3),    czone(4:5)

      IF( lwriter ) CLOSE(numtime) 
      !
   END SUBROUTINE timing_finalize
   

   SUBROUTINE wcurrent_info
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE wcurrent_info ***
      !! ** Purpose :  compute and write timing output file
      !!----------------------------------------------------------------------
      LOGICAL :: ll_ord
      CHARACTER(len=2048) :: clfmt            
   
      ! reorder the current list by elapse time      
      s_wrk => NULL()
      s_timer => s_timer_root
      DO
         ll_ord = .TRUE.
         s_timer => s_timer_root
         DO WHILE ( ASSOCIATED( s_timer%next ) )
         IF (.NOT. ASSOCIATED(s_timer%next)) EXIT
            IF ( s_timer%tsum_clock < s_timer%next%tsum_clock ) THEN 
               ALLOCATE(s_wrk)
               s_wrk = s_timer%next
               CALL insert  (s_timer, s_timer_root, s_wrk)
               CALL suppress(s_timer%next)            
               ll_ord = .FALSE.
               CYCLE            
            ENDIF           
         IF( ASSOCIATED(s_timer%next) ) s_timer => s_timer%next
         END DO         
         IF( ll_ord ) EXIT
      END DO
            
      ! write current info
      WRITE(numtime,*) 'Detailed timing for proc :', narea-1
      WRITE(numtime,*) '--------------------------'
      WRITE(numtime,*) 'Section             ',            &
      &   'Elapsed Time (s)  ','Elapsed Time (%)  ',   &
      &   'CPU Time(s)  ','CPU Time (%)  ','CPU/Elapsed  ','Frequency' 
      s_timer => s_timer_root  
      clfmt = '(1x,a,4x,f12.3,6x,f12.3,x,f12.3,2x,f12.3,6x,f7.3,2x,i9)'
      DO WHILE ( ASSOCIATED(s_timer) )
         WRITE(numtime,TRIM(clfmt))   s_timer%cname,   &
         &   s_timer%tsum_clock,s_timer%tsum_clock*100./t_elaps(2),            &
         &   s_timer%tsum_cpu  ,s_timer%tsum_cpu*100./t_cpu(2)    ,            &
         &   s_timer%tsum_cpu/s_timer%tsum_clock, s_timer%niter
         s_timer => s_timer%next
      END DO
      WRITE(numtime,*)
      !                  
   END SUBROUTINE wcurrent_info

#if defined key_mpp_mpi     
   SUBROUTINE waver_info
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE wcurrent_info ***
      !! ** Purpose :  compute and write averaged timing informations
      !!----------------------------------------------------------------------
      TYPE(alltimer), POINTER :: sl_timer_glob_root => NULL()
      TYPE(alltimer), POINTER :: sl_timer_glob      => NULL()
      TYPE(timer), POINTER :: sl_timer_ave_root => NULL()
      TYPE(timer), POINTER :: sl_timer_ave      => NULL()
      INTEGER :: icode
      INTEGER :: ierr
      LOGICAL :: ll_ord           
      CHARACTER(len=200) :: clfmt              
                 
      ! Initialised the global strucutre   
      ALLOCATE(sl_timer_glob_root, Stat=ierr)
      IF(ierr /= 0)THEN
         WRITE(numtime,*) 'Failed to allocate global timing structure in waver_info'
         RETURN
      END IF

      ALLOCATE(sl_timer_glob_root%cname     (jpnij), &
               sl_timer_glob_root%tsum_cpu  (jpnij), &
               sl_timer_glob_root%tsum_clock(jpnij), &
               sl_timer_glob_root%niter     (jpnij), Stat=ierr)
      IF(ierr /= 0)THEN
         WRITE(numtime,*) 'Failed to allocate global timing structure in waver_info'
         RETURN
      END IF
      sl_timer_glob_root%cname(:)       = ''
      sl_timer_glob_root%tsum_cpu(:)   = 0._wp
      sl_timer_glob_root%tsum_clock(:) = 0._wp
      sl_timer_glob_root%niter(:)      = 0
      sl_timer_glob_root%next => NULL()
      sl_timer_glob_root%prev => NULL()
      !ARPDBG - don't need to allocate a pointer that's immediately then
      !         set to point to some other object.
      !ALLOCATE(sl_timer_glob)
      !ALLOCATE(sl_timer_glob%cname     (jpnij))
      !ALLOCATE(sl_timer_glob%tsum_cpu  (jpnij))
      !ALLOCATE(sl_timer_glob%tsum_clock(jpnij))
      !ALLOCATE(sl_timer_glob%niter     (jpnij))
      sl_timer_glob => sl_timer_glob_root
      !
      IF( narea .EQ. 1 ) THEN
         ALLOCATE(sl_timer_ave_root)
         sl_timer_ave_root%cname       = ''
         sl_timer_ave_root%t_cpu      = 0._wp
         sl_timer_ave_root%t_clock    = 0._wp
         sl_timer_ave_root%tsum_cpu   = 0._wp
         sl_timer_ave_root%tsum_clock = 0._wp
         sl_timer_ave_root%tmax_cpu   = 0._wp
         sl_timer_ave_root%tmax_clock = 0._wp
         sl_timer_ave_root%tmin_cpu   = 0._wp
         sl_timer_ave_root%tmin_clock = 0._wp
         sl_timer_ave_root%tsub_cpu   = 0._wp
         sl_timer_ave_root%tsub_clock = 0._wp
         sl_timer_ave_root%ncount      = 0
         sl_timer_ave_root%ncount_rate = 0
         sl_timer_ave_root%ncount_max  = 0
         sl_timer_ave_root%niter       = 0
         sl_timer_ave_root%l_tdone  = .FALSE.
         sl_timer_ave_root%next => NULL()
         sl_timer_ave_root%prev => NULL()
         ALLOCATE(sl_timer_ave)
         sl_timer_ave => sl_timer_ave_root            
      ENDIF 

      ! Gather info from all processors
      s_timer => s_timer_root
      DO WHILE ( ASSOCIATED(s_timer) )
         CALL MPI_GATHER(s_timer%cname     , 20, MPI_CHARACTER,   &
                         sl_timer_glob%cname, 20, MPI_CHARACTER,   &
                         0, MPI_COMM_OPA, icode)
         CALL MPI_GATHER(s_timer%tsum_clock     , 1, MPI_DOUBLE_PRECISION,   &
                         sl_timer_glob%tsum_clock, 1, MPI_DOUBLE_PRECISION,   &
                         0, MPI_COMM_OPA, icode)
         CALL MPI_GATHER(s_timer%tsum_cpu     , 1, MPI_DOUBLE_PRECISION,   &
                         sl_timer_glob%tsum_cpu, 1, MPI_DOUBLE_PRECISION,   &
                         0, MPI_COMM_OPA, icode)
         CALL MPI_GATHER(s_timer%niter     , 1, MPI_INTEGER,   &
                         sl_timer_glob%niter, 1, MPI_INTEGER,   &
                         0, MPI_COMM_OPA, icode)

         IF( narea == 1 .AND. ASSOCIATED(s_timer%next) ) THEN
            ALLOCATE(sl_timer_glob%next)
            ALLOCATE(sl_timer_glob%next%cname     (jpnij))
            ALLOCATE(sl_timer_glob%next%tsum_cpu  (jpnij))
            ALLOCATE(sl_timer_glob%next%tsum_clock(jpnij))
            ALLOCATE(sl_timer_glob%next%niter     (jpnij))
            sl_timer_glob%next%prev => sl_timer_glob
            sl_timer_glob%next%next => NULL()
            sl_timer_glob           => sl_timer_glob%next
         ENDIF              
         s_timer => s_timer%next
      END DO      

      IF( narea == 1 ) THEN    
         ! Compute some stats
         sl_timer_glob => sl_timer_glob_root
         DO WHILE( ASSOCIATED(sl_timer_glob) )
            sl_timer_ave%cname  = sl_timer_glob%cname(1)
            sl_timer_ave%tsum_cpu   = SUM   (sl_timer_glob%tsum_cpu  (:)) / jpnij
            sl_timer_ave%tsum_clock = SUM   (sl_timer_glob%tsum_clock(:)) / jpnij
            sl_timer_ave%tmax_cpu   = MAXVAL(sl_timer_glob%tsum_cpu  (:))
            sl_timer_ave%tmax_clock = MAXVAL(sl_timer_glob%tsum_clock(:))
            sl_timer_ave%tmin_cpu   = MINVAL(sl_timer_glob%tsum_cpu  (:))
            sl_timer_ave%tmin_clock = MINVAL(sl_timer_glob%tsum_clock(:))
            sl_timer_ave%niter      = SUM   (sl_timer_glob%niter     (:))
            !
            IF( ASSOCIATED(sl_timer_glob%next) ) THEN
               ALLOCATE(sl_timer_ave%next)          
               sl_timer_ave%next%prev => sl_timer_ave
               sl_timer_ave%next%next => NULL()           
               sl_timer_ave           => sl_timer_ave%next
            ENDIF
            sl_timer_glob => sl_timer_glob%next                                
         END DO

         ! reorder the averaged list by CPU time      
         s_wrk => NULL()
         sl_timer_ave => sl_timer_ave_root
         DO
            ll_ord = .TRUE.
            sl_timer_ave => sl_timer_ave_root
            DO WHILE( ASSOCIATED( sl_timer_ave%next ) )

               IF( .NOT. ASSOCIATED(sl_timer_ave%next) ) EXIT

               IF ( sl_timer_ave%tsum_clock < sl_timer_ave%next%tsum_clock ) THEN 
                  ALLOCATE(s_wrk)
                  ! Copy data into the new object pointed to by s_wrk
                  s_wrk = sl_timer_ave%next
                  ! Insert this new timer object before our current position
                  CALL insert  (sl_timer_ave, sl_timer_ave_root, s_wrk)
                  ! Remove the old object from the list
                  CALL suppress(sl_timer_ave%next)            
                  ll_ord = .FALSE.
                  CYCLE            
               ENDIF           
               IF( ASSOCIATED(sl_timer_ave%next) ) sl_timer_ave => sl_timer_ave%next
            END DO         
            IF( ll_ord ) EXIT
         END DO

         ! write averaged info
         WRITE(numtime,"('Averaged timing on all processors :')")
         WRITE(numtime,"('-----------------------------------')")
         WRITE(numtime,"('Section',13x,'Elap. Time(s)',2x,'Elap. Time(%)',2x, &
         &   'CPU Time(s)',2x,'CPU Time(%)',2x,'CPU/Elap',1x,   &
         &   'Max elap(%)',2x,'Min elap(%)',2x,            &           
         &   'Freq')")
         sl_timer_ave => sl_timer_ave_root  
         clfmt = '((A),E15.7,2x,f6.2,5x,f12.2,5x,f6.2,5x,f7.2,2x,f12.2,4x,f6.2,2x,f9.2)'
         DO WHILE ( ASSOCIATED(sl_timer_ave) )
            WRITE(numtime,TRIM(clfmt))   sl_timer_ave%cname(1:18),                            &
            &   sl_timer_ave%tsum_clock,sl_timer_ave%tsum_clock*100.*jpnij/tot_etime,   &
            &   sl_timer_ave%tsum_cpu  ,sl_timer_ave%tsum_cpu*100.*jpnij/tot_ctime  ,   &
            &   sl_timer_ave%tsum_cpu/sl_timer_ave%tsum_clock,                          &
            &   sl_timer_ave%tmax_clock*100.*jpnij/tot_etime,                           &
            &   sl_timer_ave%tmin_clock*100.*jpnij/tot_etime,                           &                                               
            &   sl_timer_ave%niter/REAL(jpnij)
            sl_timer_ave => sl_timer_ave%next
         END DO
         WRITE(numtime,*)
         !
         DEALLOCATE(sl_timer_ave_root)
      ENDIF
      !
      DEALLOCATE(sl_timer_glob_root)
      !                  
   END SUBROUTINE waver_info
  
  
   SUBROUTINE wmpi_info
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE wmpi_time  ***
      !! ** Purpose :   compute and write a summary of MPI infos 
      !!----------------------------------------------------------------------   
      !   
      INTEGER                            :: idum, icode
      INTEGER, ALLOCATABLE, DIMENSION(:) :: iall_rank
      REAL(wp) :: ztot_ratio
      REAL(wp) :: zmax_etime, zmax_ctime, zmax_ratio, zmin_etime, zmin_ctime, zmin_ratio
      REAL(wp) :: zavg_etime, zavg_ctime, zavg_ratio
      REAL(wp), ALLOCATABLE, DIMENSION(:) :: zall_ratio
      CHARACTER(LEN=128), dimension(8) :: cllignes
      CHARACTER(LEN=128)               :: clhline, clstart_date, clfinal_date
      CHARACTER(LEN=2048)              :: clfmt    
   
      ! Gather all times
      ALLOCATE( zall_ratio(jpnij), iall_rank(jpnij) )
      IF( narea == 1 ) THEN
         iall_rank(:) = (/ (idum,idum=0,jpnij-1) /)
   
         ! Compute elapse user time
         zavg_etime = tot_etime/REAL(jpnij,wp)
         zmax_etime = MAXVAL(all_etime(:))
         zmin_etime = MINVAL(all_etime(:))

         ! Compute CPU user time
         zavg_ctime = tot_ctime/REAL(jpnij,wp)
         zmax_ctime = MAXVAL(all_ctime(:))
         zmin_ctime = MINVAL(all_ctime(:))
   
         ! Compute cpu/elapsed ratio
         zall_ratio(:) = all_ctime(:) / all_etime(:)
         ztot_ratio    = SUM(zall_ratio(:))
         zavg_ratio    = ztot_ratio/REAL(jpnij,wp)
         zmax_ratio    = MAXVAL(zall_ratio(:))
         zmin_ratio    = MINVAL(zall_ratio(:))   
   
         ! Output Format
         clhline    ='1x,13("-"),"|",18("-"),"|",14("-"),"|",18("-"),/,'
         cllignes(1)='(1x,"MPI summary report :",/,'
         cllignes(2)='1x,"--------------------",//,'
         cllignes(3)='1x,"Process Rank |"," Elapsed Time (s) |"," CPU Time (s) |"," Ratio CPU/Elapsed",/,'
         cllignes(4)='    (1x,i4,9x,"|",f12.3,6x,"|",f12.3,2x,"|",4x,f7.3,/),'
         WRITE(cllignes(4)(1:5),'(I5)') jpnij
         cllignes(5)='1x,"Total        |",f12.3,6x,"|",F12.3,2x,"|",4x,f7.3,/,'
         cllignes(6)='1x,"Minimum      |",f12.3,6x,"|",F12.3,2x,"|",4x,f7.3,/,'
         cllignes(7)='1x,"Maximum      |",f12.3,6x,"|",F12.3,2x,"|",4x,f7.3,/,'
         cllignes(8)='1x,"Average      |",f12.3,6x,"|",F12.3,2x,"|",4x,f7.3)'
         clfmt=TRIM(cllignes(1))// TRIM(cllignes(2))//TRIM(cllignes(3))//          &
           & TRIM(clhline)//TRIM(cllignes(4))//TRIM(clhline)//TRIM(cllignes(5))//  &
           & TRIM(clhline)//TRIM(cllignes(6))//TRIM(clhline)//TRIM(cllignes(7))//  &
           & TRIM(clhline)//TRIM(cllignes(8))
         WRITE(numtime, TRIM(clfmt)) &
             (iall_rank(idum),all_etime(idum),all_ctime(idum),zall_ratio(idum),idum=1, jpnij), &
             tot_etime,     tot_ctime,     ztot_ratio,   &
             zmin_etime,    zmin_ctime,    zmin_ratio,   &
             zmax_etime,    zmax_ctime,    zmax_ratio,   &
             zavg_etime,    zavg_ctime,    zavg_ratio
         WRITE(numtime,*)    
      END IF
      !
      DEALLOCATE(zall_ratio, iall_rank)
      !
   END SUBROUTINE wmpi_info
#endif   


   SUBROUTINE timing_ini_var(cdinfo)
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE timing_ini_var  ***
      !! ** Purpose :   create timing structure 
      !!----------------------------------------------------------------------
      CHARACTER(len=*), INTENT(in) :: cdinfo
      LOGICAL :: ll_section
       
      !
      IF( .NOT. ASSOCIATED(s_timer_root) ) THEN
         ALLOCATE(s_timer_root)
         s_timer_root%cname       = cdinfo
         s_timer_root%t_cpu      = 0._wp
         s_timer_root%t_clock    = 0._wp
         s_timer_root%tsum_cpu   = 0._wp
         s_timer_root%tsum_clock = 0._wp
         s_timer_root%tmax_cpu   = 0._wp
         s_timer_root%tmax_clock = 0._wp
         s_timer_root%tmin_cpu   = 0._wp
         s_timer_root%tmin_clock = 0._wp
         s_timer_root%tsub_cpu   = 0._wp
         s_timer_root%tsub_clock = 0._wp
         s_timer_root%ncount      = 0
         s_timer_root%ncount_rate = 0
         s_timer_root%ncount_max  = 0
         s_timer_root%niter       = 0
         s_timer_root%l_tdone  = .FALSE.
         s_timer_root%next => NULL()
         s_timer_root%prev => NULL()
         s_timer => s_timer_root
         !
         ALLOCATE(s_wrk)
         s_wrk => NULL()
         
      ELSE
         s_timer => s_timer_root
         ! case of already existing area (typically inside a loop)
         DO WHILE( ASSOCIATED(s_timer) ) 
            IF( TRIM(s_timer%cname) .EQ. TRIM(cdinfo) ) RETURN
            s_timer => s_timer%next
         END DO
         
         ! end of the chain
         s_timer => s_timer_root
         DO WHILE( ASSOCIATED(s_timer%next) )
            s_timer => s_timer%next
         END DO
          
         ALLOCATE(s_timer%next)      
         s_timer%next%cname       = cdinfo
         s_timer%next%t_cpu      = 0._wp
         s_timer%next%t_clock    = 0._wp
         s_timer%next%tsum_cpu   = 0._wp
         s_timer%next%tsum_clock = 0._wp  
         s_timer%next%tmax_cpu   = 0._wp
         s_timer%next%tmax_clock = 0._wp
         s_timer%next%tmin_cpu   = 0._wp
         s_timer%next%tmin_clock = 0._wp
         s_timer%next%tsub_cpu   = 0._wp
         s_timer%next%tsub_clock = 0._wp
         s_timer%next%ncount      = 0
         s_timer%next%ncount_rate = 0
         s_timer%next%ncount_max  = 0
         s_timer%next%niter       = 0
         s_timer%next%l_tdone  = .FALSE.
         s_timer%next%parent_section => NULL()
         s_timer%next%prev => s_timer
         s_timer%next%next => NULL()
         s_timer => s_timer%next

         ! are we inside a section
         s_wrk => s_timer%prev
         ll_section = .FALSE.
         DO WHILE( ASSOCIATED(s_wrk) .AND. .NOT. ll_section )
            IF( .NOT. s_wrk%l_tdone ) THEN
               ll_section = .TRUE.
               s_timer%parent_section => s_wrk 
            ENDIF
            s_wrk => s_wrk%prev
         END DO 
      ENDIF         
      !
   END SUBROUTINE timing_ini_var


   SUBROUTINE timing_reset
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE timing_reset  ***
      !! ** Purpose :   go to root of timing tree 
      !!----------------------------------------------------------------------
      l_initdone = .TRUE. 
!      IF(lwp) WRITE(numout,*)
!      IF(lwp) WRITE(numout,*) 'timing_reset : instrumented routines for timing'
!      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~'
!      CALL timing_list(s_timer_root)
!      WRITE(numout,*)
      !
   END SUBROUTINE timing_reset


   RECURSIVE SUBROUTINE timing_list(ptr)
   
      TYPE(timer), POINTER, INTENT(inout) :: ptr
      !
      IF( ASSOCIATED(ptr%next) ) CALL timing_list(ptr%next)
      IF(lwp) WRITE(numout,*)'   ', ptr%cname   
      !
   END SUBROUTINE timing_list


   SUBROUTINE insert(sd_current, sd_root ,sd_ptr)
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE insert  ***
      !! ** Purpose :   insert an element in timer structure
      !!----------------------------------------------------------------------
      TYPE(timer), POINTER, INTENT(inout) :: sd_current, sd_root, sd_ptr
      !
     
      IF( ASSOCIATED( sd_current, sd_root ) ) THEN
         ! If our current element is the root element then
         ! replace it with the one being inserted
         sd_root => sd_ptr
      ELSE
         sd_current%prev%next => sd_ptr
      END IF
      sd_ptr%next     => sd_current
      sd_ptr%prev     => sd_current%prev
      sd_current%prev => sd_ptr
      ! Nullify the pointer to the new element now that it is held
      ! within the list. If we don't do this then a subsequent call
      ! to ALLOCATE memory to this pointer will fail.
      sd_ptr => NULL()
      !    
   END SUBROUTINE insert
  
  
   SUBROUTINE suppress(sd_ptr)
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE suppress  ***
      !! ** Purpose :   supress an element in timer structure
      !!----------------------------------------------------------------------
      TYPE(timer), POINTER, INTENT(inout) :: sd_ptr
      !
      TYPE(timer), POINTER :: sl_temp
    
      sl_temp => sd_ptr
      sd_ptr => sd_ptr%next    
      IF ( ASSOCIATED(sl_temp%next) ) sl_temp%next%prev => sl_temp%prev
      DEALLOCATE(sl_temp)
      sl_temp => NULL()
      !
    END SUBROUTINE suppress

   !!=====================================================================
END MODULE timing
