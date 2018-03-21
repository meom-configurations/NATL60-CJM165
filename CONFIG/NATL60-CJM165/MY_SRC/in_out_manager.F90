MODULE in_out_manager   
   !!======================================================================
   !!                       ***  MODULE  in_out_manager  ***
   !! I/O manager utilities : Defines run parameters together with logical units
   !!=====================================================================
   !! History :   1.0  !  2002-06  (G. Madec)   original code
   !!             2.0  !  2006-07  (S. Masson)  iom, add ctl_stop, ctl_warn
   !!             3.0  !  2008-06  (G. Madec)   add ctmp4 to ctmp10
   !!             3.2  !  2009-08  (S. MAsson)  add new ctl_opn
   !!             3.3  !  2010-10  (A. Coward)  add NetCDF4 usage
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   USE par_oce       ! ocean parameter
   USE lib_print     ! formated print library
   USE nc4interface  ! NetCDF4 interface

   IMPLICIT NONE
   PUBLIC

 
   !
   !!----------------------------------------------------------------------
   !!                   namrun namelist parameters
   !!----------------------------------------------------------------------
   CHARACTER(lc) ::   cn_exp           !: experiment name used for output filename
   CHARACTER(lc) ::   cn_ocerst_in     !: suffix of ocean restart name (input)
   CHARACTER(lc) ::   cn_ocerst_indir  !: restart input directory
   CHARACTER(lc) ::   cn_ocerst_out    !: suffix of ocean restart name (output)
   CHARACTER(lc) ::   cn_ocerst_outdir !: restart output directory
   CHARACTER(len=255) ::  cn_dirout = "./"        !: prefix for model outputfile
   LOGICAL       ::   ln_rstart        !: start from (F) rest or (T) a restart file
   LOGICAL       ::   ln_rst_list      !: output restarts at list of times (T) or by frequency (F)
   INTEGER       ::   nn_no            !: job number
   INTEGER       ::   nn_rstctl        !: control of the time step (0, 1 or 2)
   INTEGER       ::   nn_rstssh   = 0  !: hand made initilization of ssh or not (1/0)
   INTEGER       ::   nn_it000         !: index of the first time step
   INTEGER       ::   nn_itend         !: index of the last time step
   INTEGER       ::   nn_date0         !: initial calendar date aammjj
   INTEGER       ::   nn_leapy         !: Leap year calendar flag (0/1 or 30)
   INTEGER       ::   nn_istate        !: initial state output flag (0/1)
   INTEGER       ::   nn_write         !: model standard output frequency
   INTEGER       ::   nn_stock         !: restart file frequency
   INTEGER, DIMENSION(10) :: nn_stocklist  !: restart dump times
   LOGICAL       ::   ln_dimgnnn       !: type of dimgout. (F): 1 file for all proc
                                                       !:                  (T): 1 file per proc
   LOGICAL       ::   ln_mskland       !: mask land points in NetCDF outputs (costly: + ~15%)
   LOGICAL       ::   ln_cfmeta        !: output additional data to netCDF files required for compliance with the CF metadata standard
   LOGICAL       ::   ln_clobber       !: clobber (overwrite) an existing file
   INTEGER       ::   nn_chunksz       !: chunksize (bytes) for NetCDF file (works only with iom_nf90 routines)
#if defined key_netcdf4
   !!----------------------------------------------------------------------
   !!                   namnc4 namelist parameters                         (key_netcdf4)
   !!----------------------------------------------------------------------
   ! The following four values determine the partitioning of the output fields
   ! into netcdf4 chunks. They are unrelated to the nn_chunk_sz setting which is
   ! for runtime optimisation. The individual netcdf4 chunks can be optionally 
   ! gzipped (recommended) leading to significant reductions in I/O volumes 
   !                         !!!**  variables only used with iom_nf90 routines and key_netcdf4 **
   INTEGER ::   nn_nchunks_i   !: number of chunks required in the i-dimension 
   INTEGER ::   nn_nchunks_j   !: number of chunks required in the j-dimension 
   INTEGER ::   nn_nchunks_k   !: number of chunks required in the k-dimension 
   INTEGER ::   nn_nchunks_t   !: number of chunks required in the t-dimension 
   LOGICAL ::   ln_nc4zip      !: netcdf4 usage: (T) chunk and compress output using the HDF5 sublayers of netcdf4
   !                           !                 (F) ignore chunking request and use the netcdf4 library 
   !                           !                     to produce netcdf3-compatible files 
#endif
!$AGRIF_DO_NOT_TREAT
   TYPE(snc4_ctl)     :: snc4set        !: netcdf4 chunking control structure (always needed for decision making)
!$AGRIF_END_DO_NOT_TREAT


   !! conversion of DOCTOR norm namelist name into model name
   !! (this should disappear in a near futur)

   CHARACTER(lc) ::   cexper                      !: experiment name used for output filename
   INTEGER       ::   no                          !: job number
   INTEGER       ::   nrstdt                      !: control of the time step (0, 1 or 2)
   INTEGER       ::   nit000                      !: index of the first time step
   INTEGER       ::   nitend                      !: index of the last time step
   INTEGER       ::   ndate0                      !: initial calendar date aammjj
   INTEGER       ::   nleapy                      !: Leap year calendar flag (0/1 or 30)
   INTEGER       ::   ninist                      !: initial state output flag (0/1)
   INTEGER       ::   nwrite                      !: model standard output frequency
   INTEGER       ::   nstock                      !: restart file frequency
   INTEGER, DIMENSION(10) :: nstocklist           !: restart dump times

   !!----------------------------------------------------------------------
   !! was in restart but moved here because of the OFF line... better solution should be found...
   !!----------------------------------------------------------------------
   INTEGER ::   nitrst                !: time step at which restart file should be written
   LOGICAL ::   lrst_oce              !: logical to control the oce restart write 
   INTEGER ::   numror = 0            !: logical unit for ocean restart (read). Init to 0 is needed for SAS (in daymod.F90)
   INTEGER ::   numrow                !: logical unit for ocean restart (write)
   INTEGER ::   nrst_lst              !: number of restart to output next

   !!----------------------------------------------------------------------
   !!                    output monitoring
   !!----------------------------------------------------------------------
   LOGICAL ::   ln_ctl       !: run control for debugging
   INTEGER ::   nn_timing    !: run control for timing
   INTEGER ::   nn_print     !: level of print (0 no print)
   INTEGER ::   nn_ictls     !: Start i indice for the SUM control
   INTEGER ::   nn_ictle     !: End   i indice for the SUM control
   INTEGER ::   nn_jctls     !: Start j indice for the SUM control
   INTEGER ::   nn_jctle     !: End   j indice for the SUM control
   INTEGER ::   nn_isplt     !: number of processors following i
   INTEGER ::   nn_jsplt     !: number of processors following j
   INTEGER ::   nn_bench     !: benchmark parameter (0/1)
   INTEGER ::   nn_bit_cmp   =    0    !: bit reproducibility  (0/1)

   !                                          
   INTEGER ::   nprint, nictls, nictle, njctls, njctle, isplt, jsplt, nbench    !: OLD namelist names

   INTEGER ::   ijsplt     =    1      !: nb of local domain = nb of processors

   !!----------------------------------------------------------------------
   !!                        logical units
   !!----------------------------------------------------------------------
   INTEGER ::   numstp          =   -1      !: logical unit for time step
   INTEGER ::   numtime         =   -1      !: logical unit for timing
   INTEGER ::   numout          =    6      !: logical unit for output print; Set to stdout to ensure any early
                                            !  output can be collected; do not change
   INTEGER ::   numnam_ref      =   -1      !: logical unit for reference namelist
   INTEGER ::   numnam_cfg      =   -1      !: logical unit for configuration specific namelist
   INTEGER ::   numond          =   -1      !: logical unit for Output Namelist Dynamics
   INTEGER ::   numnam_ice_ref  =   -1      !: logical unit for ice reference namelist
   INTEGER ::   numnam_ice_cfg  =   -1      !: logical unit for ice reference namelist
   INTEGER ::   numoni          =   -1      !: logical unit for Output Namelist Ice
   INTEGER ::   numevo_ice      =   -1      !: logical unit for ice variables (temp. evolution)
   INTEGER ::   numsol          =   -1      !: logical unit for solver statistics
   INTEGER ::   numdct_in       =   -1      !: logical unit for transports computing
   INTEGER ::   numdct_vol      =   -1      !: logical unit for voulume transports output
   INTEGER ::   numdct_heat     =   -1      !: logical unit for heat    transports output
   INTEGER ::   numdct_salt     =   -1      !: logical unit for salt    transports output
   INTEGER ::   numfl           =   -1      !: logical unit for floats ascii output
   INTEGER ::   numflo          =   -1      !: logical unit for floats ascii output

   !!----------------------------------------------------------------------
   !!                          Run control  
   !!----------------------------------------------------------------------
   INTEGER       ::   nstop = 0             !: error flag (=number of reason for a premature stop run)
   INTEGER       ::   nwarn = 0             !: warning flag (=number of warning found during the run)
   CHARACTER(lc) ::   ctmp1, ctmp2, ctmp3   !: temporary characters 1 to 3
   CHARACTER(lc) ::   ctmp4, ctmp5, ctmp6   !: temporary characters 4 to 6
   CHARACTER(lc) ::   ctmp7, ctmp8, ctmp9   !: temporary characters 7 to 9
   CHARACTER(lc) ::   ctmp10                !: temporary character 10
   CHARACTER(lc) ::   cform_err = "(/,' ===>>> : E R R O R',     /,'         ===========',/)"       !:
   CHARACTER(lc) ::   cform_war = "(/,' ===>>> : W A R N I N G', /,'         ===============',/)"   !:
   LOGICAL       ::   lwm      = .FALSE.    !: boolean : true on the 1st processor only (always)
   LOGICAL       ::   lwp      = .FALSE.    !: boolean : true on the 1st processor only .OR. ln_ctl
   LOGICAL       ::   lsp_area = .TRUE.     !: to make a control print over a specific area
   CHARACTER(lc) ::   cxios_context         !: context name used in xios

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: in_out_manager.F90 5518 2015-06-30 13:11:42Z clem $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!=====================================================================
END MODULE in_out_manager
