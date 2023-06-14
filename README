#----------------------------------------------------------------------
#
#  GC2CMAQ v1.0
#  Developed on 2022/01/10 by Yiming Liu
#  Assistant Professor
#  School of Atmospheric Sciences, Sun Yat-sen University
#
#  Comments and suggestions are welcome.
#  Please mail to liuym88@mail.sysu.edu.cn
#
#  Update on 2022/10/01:  Keep the interpolated concentration above zero.
#
#
#----------------------------------------------------------------------

This programe is developed to provide IC/BC of the CMAQ model from 
the GEOS-Chem outputs. It is developed based on MOZART2CAMX 
(https://www.camx.com). Thanks.

GC2CMAQ reads one GEOS-Chem (version 13.3) output files (in GMT) and 
horizontally and vertically interpolates the data onto a CMAQ (v5.2.1)
domain to generate date-specific initial condition and/or boundary 
condition files.  

If generating CMAQ outputs, MCIP METCRO3D and METCRO2D files
are needed to determine the layer structure.  Note that only the first
hour of met data will be used to determine the layer structure.   

Horizontal interpolation is bilinear interpolation.

Vertical interpolation uses pressure levels of each layer interface.  A simple
interpolation is performed such that concentrations in each GEOS-Chem layer are
weighted based on the range in pressure within each CMAQ layer.  For
simplicity, differences in density between the adjacent GEOS-Chem layers are
assumed to be the same for the vertical weightings.

Speciation profiles are included to map the GEOS-Chem species to 
CMAQ (CB05, SAPRC07TIC).  If recompiling is
needed, select the appropriate speciation mapping file in src/G2Lconv_*.EXT and
see the top of Makefile on how to compile.

GC2CMAQ reads in 1 GEOS-Chem output files (boundary condition file) and 
outputs BCs every 3 hours (9 time periods for CMAQ). The user can specify the hour
of the initial conditions; values will be extracted at either the specified
time period or from the nearest time period BEFOREHAND.  IC outputs will have
a time stamp for the entire date, so no time shifting should be needed.
 
*****

Steps to process:
1.  Run the GEOS-Chem model and get the boundary condition output file.
    1 day per file, including 9 time periods (should attach the first hour data
    from tomorrow's file to today's file). 
    Please contact Yiming Liu (liuym88@mail.sysu.edu.cn) or 
    Xiao Lu (luxiao25@mail.sysu.edu.cn) from Sun Yat-sen University for the 
    available GEOS-Chem output data. The data can be accessed upon request.

2.  Build gc2cmaq executable (src folder). 
    Select the CMAQ mechanism and compile it (make SAPRC07TIC_AE6I).

3.  Run gc2cmaq (script/run.gc2cmaq) to generate CMAQ IC/BCs.
    Outputs are in GMT.
________________________________________________________________________________

Sample job script for CMAQ:

setenv OUTFILEBC  $OUTPATH/bc.${GRID_NAME}.cmaq.saprc07tic_ae6i_aqkmti.${YYYYMMDD}.ncf
setenv OUTFILEIC  $OUTPATH/ic.${GRID_NAME}.cmaq.saprc07tic_ae6i_aqkmti.${YYYYMMDD}.hr0.ncf
setenv INFILE1    $DPATH/input/GEOSChem.BoundaryConditions.${YYYYMMDD}.nc4
setenv INFILEMET3D  $DPATH/mcip/${GRID_NAME}/2020234/METCRO3D_2020234
setenv INFILEMET2D  $DPATH/mcip/${GRID_NAME}/2020234/METCRO2D_2020234

../src/gc2cmaq_SAPRC07TIC_AE6 << IEOF 
ProcessDateYYYYMMDD|$YYYYMMDD
Output BC file?    |.true.
Output IC file?    |.true.
If IC, starting hr |0
IEOF

line 1: Gregorian date to process (YYYYMMDD) 
line 2: Option to output a boundary condition file (true/false)
line 3: Option to output an initial conditions file (true/false)
line 4: If outputting an IC file, hour to extract (GMT)
        Note: GEOS-Chem outputs are 3 hourly so the closest time BEFORE the 
        specified hour will be used if the times do not match

Job script only needs first 4 lines since the rest can be pulled
from the MCIP files
