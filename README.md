# NATL60-CJM165
This repository holds the source code corresponding to NEMO configuration for NATL60-CJM165.
This configuration is based on rev_6355 of NEMO_v3_STABLE branch and was operated with XIOS rev 703 of XIOS/branchs/xios-1.0.  NEMO and XIOS can be downloaded from the IPSL forge with the commands:

svn co -r 6355 http://forge.ipsl.jussieu.fr/nemo/svn/branches/2015/nemo_v3_6_STABLE/NEMOGCM

svn co -r 703 http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/branchs/xios-1.0

## Organisation of the repository
In the CONFIG directory we give the NATL60-CJM165 holding the needed fcm files (cpp and arch) and  classical  MY_SRC/ and EXP00/ sub-directories.

  MY_SRC holds fortran modules modified with respect to the reference NEMO code
  
  EXP00 holds namelist for the ocean an ice, as well as the xml files suitable for XIOS output.
