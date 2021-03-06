# NATL60-CJM165
This repository holds the source code corresponding to NEMO configuration for NATL60-CJM165.
This configuration is based on rev_6355 of NEMO_v3_STABLE branch and was operated with XIOS rev_703 of XIOS/branchs/xios-1.0.  NEMO and XIOS can be downloaded from the IPSL forge with the commands:

```svn co -r 6355 http://forge.ipsl.jussieu.fr/nemo/svn/branches/2015/nemo_v3_6_STABLE/NEMOGCM```

```svn co -r 703 http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/branchs/xios-1.0```

## Description of the repository:
  This repository holds all the information needed to build the numerical code used for NATL60-CJM165 configuration, from the reference versions indicated above.
  
  It provides the ```CONFIG/NATL60-CJM165``` directory to be used for NATL60-CJM165:
  
* ```CONFIG/NATL60-CJM165 ``` : fcm files for cpp and compilation options.
* ```CONFIG/NATL60-CJM165/MY_SRC/``` : fortran code differing from the reference NEMO code.
* ```CONFIG/NATL60-CJM165/EXP00/``` : namelist files (ocean, ice and top), xml files for xios output.

## Run-time files:
  Most of the run-time files are indicated in the namelist files, except for :
* bathymetry : ```NATL60_bathymeter_zps_gebco_v4.nc```
* coordinates : ```NATL60_coordinates_v4.nc```
* ice initialization : ```Init_Ice_GLORYS1V1_NSIDC_BOOTSTRAP_y1989m01_new.nc```

Acces to run-time files can be made available on request.
