MODULE mppini
   !!==============================================================================
   !!                       ***  MODULE mppini   ***
   !! Ocean initialization : distributed memory computing initialization
   !!==============================================================================

   !!----------------------------------------------------------------------
   !!   mpp_init       : Lay out the global domain over processors
   !!   mpp_init2      : Lay out the global domain over processors 
   !!                    with land processor elimination
   !!   mpp_init_ioispl: IOIPSL initialization in mpp
   !!----------------------------------------------------------------------
   !! * Modules used
   USE dom_oce         ! ocean space and time domain 
   USE in_out_manager  ! I/O Manager
   USE lib_mpp         ! distribued memory computing library
   USE ioipsl

   IMPLICIT NONE
   PRIVATE

   PUBLIC mpp_init       ! called by opa.F90
   PUBLIC mpp_init2      ! called by opa.F90

   !! * Substitutions
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: mppini.F90 4679 2014-06-20 10:17:06Z epico $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

#if ! defined key_mpp_mpi
   !!----------------------------------------------------------------------
   !!   Default option :                            shared memory computing
   !!----------------------------------------------------------------------

   SUBROUTINE mpp_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE mpp_init  ***
      !!
      !! ** Purpose :   Lay out the global domain over processors.
      !!
      !! ** Method  :   Shared memory computing, set the local processor
      !!      variables to the value of the global domain
      !!
      !! History :
      !!   9.0  !  04-01  (G. Madec, J.M. Molines)  F90 : free form, north fold jpni >1
      !!----------------------------------------------------------------------

      ! No mpp computation
      nimpp  = 1
      njmpp  = 1
      nlci   = jpi
      nlcj   = jpj
      nldi   = 1
      nldj   = 1
      nlei   = jpi
      nlej   = jpj
      nperio = jperio
      nbondi = 2
      nbondj = 2
      nidom  = FLIO_DOM_NONE
      npolj = jperio

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'mpp_init(2) : NO massively parallel processing'
         WRITE(numout,*) '~~~~~~~~~~~: '
         WRITE(numout,*) '         nperio = ', nperio
         WRITE(numout,*) '         npolj  = ', npolj
         WRITE(numout,*) '         nimpp  = ', nimpp
         WRITE(numout,*) '         njmpp  = ', njmpp
      ENDIF

      IF(  jpni /= 1 .OR. jpnj /= 1 .OR. jpnij /= 1 ) &
          CALL ctl_stop( 'equality  jpni = jpnj = jpnij = 1 is not satisfied',   &
          &              'the domain is lay out for distributed memory computing! ' )

   END SUBROUTINE mpp_init


   SUBROUTINE mpp_init2 
      CALL mpp_init                             ! same routine as mpp_init
   END SUBROUTINE mpp_init2

#else
   !!----------------------------------------------------------------------
   !!   'key_mpp_mpi'          OR         MPI massively parallel processing
   !!----------------------------------------------------------------------

   SUBROUTINE mpp_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE mpp_init  ***
      !!                    
      !! ** Purpose :   Lay out the global domain over processors.
      !!
      !! ** Method  :   Global domain is distributed in smaller local domains.
      !!      Periodic condition is a function of the local domain position
      !!      (global boundary or neighbouring domain) and of the global
      !!      periodic
      !!      Type :         jperio global periodic condition
      !!                     nperio local  periodic condition
      !!
      !! ** Action  : - set domain parameters
      !!                    nimpp     : longitudinal index 
      !!                    njmpp     : latitudinal  index
      !!                    nperio    : lateral condition type 
      !!                    narea     : number for local area
      !!                    nlci      : first dimension
      !!                    nlcj      : second dimension
      !!                    nbondi    : mark for "east-west local boundary"
      !!                    nbondj    : mark for "north-south local boundary"
      !!                    nproc     : number for local processor
      !!                    noea      : number for local neighboring processor
      !!                    nowe      : number for local neighboring processor
      !!                    noso      : number for local neighboring processor
      !!                    nono      : number for local neighboring processor
      !!
      !! History :
      !!        !  94-11  (M. Guyon)  Original code
      !!        !  95-04  (J. Escobar, M. Imbard)
      !!        !  98-02  (M. Guyon)  FETI method
      !!        !  98-05  (M. Imbard, J. Escobar, L. Colombet )  SHMEM and MPI versions
      !!   8.5  !  02-08  (G. Madec)  F90 : free form
      !!   3.4  !  11-11  (C. Harris) decomposition changes for running with CICE
      !!----------------------------------------------------------------------
      INTEGER  ::   ji, jj, jn   ! dummy loop indices
      INTEGER  ::   ii, ij, ifreq, il1, il2            ! local integers
      INTEGER  ::   iresti, irestj, ijm1, imil, inum   !   -      -
      REAL(wp) ::   zidom, zjdom                       ! local scalars
      INTEGER, DIMENSION(jpni,jpnj) ::   iimppt, ijmppt, ilcit, ilcjt   ! local workspace
      !!----------------------------------------------------------------------

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'mpp_init : Message Passing MPI'
      IF(lwp) WRITE(numout,*) '~~~~~~~~'


      !  1. Dimension arrays for subdomains
      ! -----------------------------------
      !  Computation of local domain sizes ilcit() ilcjt()
      !  These dimensions depend on global sizes jpni,jpnj and jpiglo,jpjglo
      !  The subdomains are squares leeser than or equal to the global
      !  dimensions divided by the number of processors minus the overlap
      !  array (cf. par_oce.F90).
      
      nreci  = 2 * jpreci
      nrecj  = 2 * jprecj
      iresti = MOD( jpiglo - nreci , jpni )
      irestj = MOD( jpjglo - nrecj , jpnj )

      IF(  iresti == 0 )   iresti = jpni

#if defined key_nemocice_decomp
      ! In order to match CICE the size of domains in NEMO has to be changed
      ! The last line of blocks (west) will have fewer points

      DO jj = 1, jpnj
         DO ji=1, jpni-1
            ilcit(ji,jj) = jpi
         END DO
         ilcit(jpni,jj) = jpiglo - (jpni - 1) * (jpi - nreci)
      END DO

#else

      DO jj = 1, jpnj
         DO ji = 1, iresti
            ilcit(ji,jj) = jpi
         END DO
         DO ji = iresti+1, jpni
            ilcit(ji,jj) = jpi -1
         END DO
      END DO
      
#endif
      nfilcit(:,:) = ilcit(:,:)
      IF( irestj == 0 )   irestj = jpnj

#if defined key_nemocice_decomp
      ! Same change to domains in North-South direction as in East-West. 
      DO ji=1,jpni
         DO jj=1,jpnj-1
            ilcjt(ji,jj) = jpj
         END DO
         ilcjt(ji,jpnj) = jpjglo - (jpnj - 1) * (jpj - nrecj)
      END DO

#else

      DO ji = 1, jpni
         DO jj = 1, irestj
            ilcjt(ji,jj) = jpj
         END DO
         DO jj = irestj+1, jpnj
            ilcjt(ji,jj) = jpj -1
         END DO
      END DO
      
#endif
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) '           defines mpp subdomains'
         WRITE(numout,*) '           ----------------------'
         WRITE(numout,*) '           iresti=',iresti,' irestj=',irestj
         WRITE(numout,*) '           jpni  =',jpni  ,' jpnj  =',jpnj
         ifreq = 4
         il1   = 1
         DO jn = 1, (jpni-1)/ifreq+1
            il2 = MIN( jpni, il1+ifreq-1 )
            WRITE(numout,*)
            WRITE(numout,9200) ('***',ji = il1,il2-1)
            DO jj = jpnj, 1, -1
               WRITE(numout,9203) ('   ',ji = il1,il2-1)
               WRITE(numout,9202) jj, ( ilcit(ji,jj),ilcjt(ji,jj),ji = il1,il2 )
               WRITE(numout,9203) ('   ',ji = il1,il2-1)
               WRITE(numout,9200) ('***',ji = il1,il2-1)
            END DO
            WRITE(numout,9201) (ji,ji = il1,il2)
            il1 = il1+ifreq
         END DO
 9200    FORMAT('     ***',20('*************',a3))
 9203    FORMAT('     *     ',20('         *   ',a3))
 9201    FORMAT('        ',20('   ',i3,'          '))
 9202    FORMAT(' ',i3,' *  ',20(i3,'  x',i3,'   *   '))
      ENDIF

      zidom = nreci
      DO ji = 1, jpni
         zidom = zidom + ilcit(ji,1) - nreci
      END DO
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*)' sum ilcit(i,1) = ', zidom, ' jpiglo = ', jpiglo
      
      zjdom = nrecj
      DO jj = 1, jpnj
         zjdom = zjdom + ilcjt(1,jj) - nrecj
      END DO
      IF(lwp) WRITE(numout,*)' sum ilcit(1,j) = ', zjdom, ' jpjglo = ', jpjglo
      IF(lwp) WRITE(numout,*)
      

      !  2. Index arrays for subdomains
      ! -------------------------------
      
      iimppt(:,:) = 1
      ijmppt(:,:) = 1
      
      IF( jpni > 1 ) THEN
         DO jj = 1, jpnj
            DO ji = 2, jpni
               iimppt(ji,jj) = iimppt(ji-1,jj) + ilcit(ji-1,jj) - nreci
            END DO
         END DO
      ENDIF
      nfiimpp(:,:)=iimppt(:,:)

      IF( jpnj > 1 ) THEN
         DO jj = 2, jpnj
            DO ji = 1, jpni
               ijmppt(ji,jj) = ijmppt(ji,jj-1)+ilcjt(ji,jj-1)-nrecj
            END DO
         END DO
      ENDIF
      
      ! 3. Subdomain description
      ! ------------------------

      DO jn = 1, jpnij
         ii = 1 + MOD( jn-1, jpni )
         ij = 1 + (jn-1) / jpni
         nfipproc(ii,ij) = jn - 1
         nimppt(jn) = iimppt(ii,ij)
         njmppt(jn) = ijmppt(ii,ij)
         nlcit (jn) = ilcit (ii,ij)     
         nlci       = nlcit (jn)     
         nlcjt (jn) = ilcjt (ii,ij)     
         nlcj       = nlcjt (jn)
         nbondj = -1                                   ! general case
         IF( jn   >  jpni          )   nbondj = 0      ! first row of processor
         IF( jn   >  (jpnj-1)*jpni )   nbondj = 1      ! last  row of processor
         IF( jpnj == 1             )   nbondj = 2      ! one processor only in j-direction
         ibonjt(jn) = nbondj
         
         nbondi = 0                                    ! 
         IF( MOD( jn, jpni ) == 1 )   nbondi = -1      !
         IF( MOD( jn, jpni ) == 0 )   nbondi =  1      !
         IF( jpni            == 1 )   nbondi =  2      ! one processor only in i-direction
         ibonit(jn) = nbondi
         
         nldi =  1   + jpreci
         nlei = nlci - jpreci
         IF( nbondi == -1 .OR. nbondi == 2 )   nldi = 1
         IF( nbondi ==  1 .OR. nbondi == 2 )   nlei = nlci
         nldj =  1   + jprecj
         nlej = nlcj - jprecj
         IF( nbondj == -1 .OR. nbondj == 2 )   nldj = 1
         IF( nbondj ==  1 .OR. nbondj == 2 )   nlej = nlcj
         nldit(jn) = nldi
         nleit(jn) = nlei
         nldjt(jn) = nldj
         nlejt(jn) = nlej
      END DO
      

      ! 4. From global to local
      ! -----------------------

      nperio = 0
      IF( jperio == 2 .AND. nbondj == -1 )   nperio = 2


      ! 5. Subdomain neighbours
      ! ----------------------

      nproc = narea - 1
      noso  = nproc - jpni
      nowe  = nproc - 1
      noea  = nproc + 1
      nono  = nproc + jpni

      nlcj = nlcjt(narea)  
      nlci = nlcit(narea)  
      nldi = nldit(narea)
      nlei = nleit(narea)
      nldj = nldjt(narea)
      nlej = nlejt(narea)
      nbondi = ibonit(narea)
      nbondj = ibonjt(narea)
      nimpp  = nimppt(narea)  
      njmpp  = njmppt(narea)  

     ! Save processor layout in layout.dat file 
       IF (lwp) THEN
        CALL ctl_opn( inum, 'layout.dat', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE., narea )
        WRITE(inum,'(a)') '   jpnij     jpi     jpj     jpk  jpiglo  jpjglo'
        WRITE(inum,'(6i8)') jpnij,jpi,jpj,jpk,jpiglo,jpjglo
        WRITE(inum,'(a)') 'RANK nlci nlcj nldi nldj nlei nlej nimpp njmpp'

        DO  jn = 1, jpnij
         WRITE(inum,'(9i5)') jn-1, nlcit(jn), nlcjt(jn), &
                                      nldit(jn), nldjt(jn), &
                                      nleit(jn), nlejt(jn), &
                                      nimppt(jn), njmppt(jn)
        END DO
        CLOSE(inum)   
      END IF


      ! w a r n i n g  narea (zone) /= nproc (processors)!

      IF( jperio == 1 .OR. jperio == 4 .OR. jperio == 6 ) THEN
         IF( jpni == 1 )THEN
            nbondi = 2
            nperio = 1
         ELSE
            nbondi = 0
         ENDIF
         IF( MOD( narea, jpni ) == 0 ) THEN
            noea = nproc-(jpni-1)
            npne = npne-jpni
            npse = npse-jpni
         ENDIF
         IF( MOD( narea, jpni ) == 1 ) THEN
            nowe = nproc+(jpni-1)
            npnw = npnw+jpni
            npsw = npsw+jpni
         ENDIF
         nbsw = 1
         nbnw = 1
         nbse = 1
         nbne = 1
         IF( nproc < jpni ) THEN
            nbsw = 0
            nbse = 0
         ENDIF
         IF( nproc >= (jpnj-1)*jpni ) THEN
            nbnw = 0
            nbne = 0
         ENDIF
      ENDIF
      npolj = 0
      IF( jperio == 3 .OR. jperio == 4 ) THEN
         ijm1 = jpni*(jpnj-1)
         imil = ijm1+(jpni+1)/2
         IF( narea > ijm1 ) npolj = 3
         IF( MOD(jpni,2) == 1 .AND. narea == imil ) npolj = 4
         IF( npolj == 3 ) nono = jpni*jpnj-narea+ijm1
      ENDIF
      IF( jperio == 5 .OR. jperio == 6 ) THEN
          ijm1 = jpni*(jpnj-1)
          imil = ijm1+(jpni+1)/2
          IF( narea > ijm1) npolj = 5
          IF( MOD(jpni,2) == 1 .AND. narea == imil ) npolj = 6
          IF( npolj == 5 ) nono = jpni*jpnj-narea+ijm1
      ENDIF

      ! Periodicity : no corner if nbondi = 2 and nperio != 1

      IF(lwp) THEN
         WRITE(numout,*) ' nproc  = ', nproc
         WRITE(numout,*) ' nowe   = ', nowe  , ' noea   =  ', noea
         WRITE(numout,*) ' nono   = ', nono  , ' noso   =  ', noso
         WRITE(numout,*) ' nbondi = ', nbondi
         WRITE(numout,*) ' nbondj = ', nbondj
         WRITE(numout,*) ' npolj  = ', npolj
         WRITE(numout,*) ' nperio = ', nperio
         WRITE(numout,*) ' nlci   = ', nlci
         WRITE(numout,*) ' nlcj   = ', nlcj
         WRITE(numout,*) ' nimpp  = ', nimpp
         WRITE(numout,*) ' njmpp  = ', njmpp
         WRITE(numout,*) ' nbse   = ', nbse  , ' npse   = ', npse
         WRITE(numout,*) ' nbsw   = ', nbsw  , ' npsw   = ', npsw
         WRITE(numout,*) ' nbne   = ', nbne  , ' npne   = ', npne
         WRITE(numout,*) ' nbnw   = ', nbnw  , ' npnw   = ', npnw
      ENDIF

      IF( nperio == 1 .AND. jpni /= 1 ) CALL ctl_stop( ' mpp_init: error on cyclicity' )

      ! Prepare mpp north fold

      IF (jperio >= 3 .AND. jperio <= 6 .AND. jpni > 1 ) THEN
         CALL mpp_ini_north
      END IF

      ! Prepare NetCDF output file (if necessary)
      CALL mpp_init_ioipsl

   END SUBROUTINE mpp_init

#  include "mppini_2.h90"

# if defined key_dimgout
   !!----------------------------------------------------------------------
   !!   'key_dimgout'                  NO use of NetCDF files
   !!----------------------------------------------------------------------
   SUBROUTINE mpp_init_ioipsl       ! Dummy routine
   END SUBROUTINE mpp_init_ioipsl  
# else
   SUBROUTINE mpp_init_ioipsl
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE mpp_init_ioipsl  ***
      !!
      !! ** Purpose :   
      !!
      !! ** Method  :   
      !!
      !! History :
      !!   9.0  !  04-03  (G. Madec )  MPP-IOIPSL 
      !!   " "  !  08-12  (A. Coward)  addition in case of jpni*jpnj < jpnij
      !!----------------------------------------------------------------------
      INTEGER, DIMENSION(2) ::   iglo, iloc, iabsf, iabsl, ihals, ihale, idid
      !!----------------------------------------------------------------------

      ! The domain is split only horizontally along i- or/and j- direction
      ! So we need at the most only 1D arrays with 2 elements.
      ! Set idompar values equivalent to the jpdom_local_noextra definition
      ! used in IOM. This works even if jpnij .ne. jpni*jpnj.
      iglo(1) = jpiglo
      iglo(2) = jpjglo
      iloc(1) = nlci
      iloc(2) = nlcj
      iabsf(1) = nimppt(narea)
      iabsf(2) = njmppt(narea)
      iabsl(:) = iabsf(:) + iloc(:) - 1
      ihals(1) = nldi - 1
      ihals(2) = nldj - 1
      ihale(1) = nlci - nlei
      ihale(2) = nlcj - nlej
      idid(1) = 1
      idid(2) = 2

      IF(lwp) THEN
          WRITE(numout,*)
          WRITE(numout,*) 'mpp_init_ioipsl :   iloc  = ', iloc (1), iloc (2)
          WRITE(numout,*) '~~~~~~~~~~~~~~~     iabsf = ', iabsf(1), iabsf(2)
          WRITE(numout,*) '                    ihals = ', ihals(1), ihals(2)
          WRITE(numout,*) '                    ihale = ', ihale(1), ihale(2)
      ENDIF
      !
      CALL flio_dom_set ( jpnij, nproc, idid, iglo, iloc, iabsf, iabsl, ihals, ihale, 'BOX', nidom)
      !
   END SUBROUTINE mpp_init_ioipsl  

# endif
#endif

   !!======================================================================
END MODULE mppini
