F90 = ifort
F90_OPTS = -O2 -mieee-fp -align dcommons -convert big_endian #-ipo -static-intel -openmp

IOAPI_INC  = /WORK/grmc_yyhong_1/share/software/ioapi_3.1 
IOAPI_LIB  = /WORK/grmc_yyhong_1/share/software/ioapi_3.1
NETCDF_LIB = /WORK/app/netcdf/4.3.2/01-CF-14/lib 
PROGRAM = gc2cmaq

LIBS =    -L$(IOAPI_LIB) -lioapi -L$(NETCDF_LIB) -lnetcdf -lnetcdff
INCLUDE = -I$(IOAPI_INC) \
          -I/WORK/app/netcdf/4.3.2/01-CF-14/include

#F90 = pgf90
#F90_OPTS = -O2 -Mdalign  -Bstatic_pgi -byteswapio

###F90 = ifort
###F90_OPTS = -O2 -align dcommons -convert big_endian

###F90 = gfortran
###F90_OPTS = -O2 -fconvert=big-endian -frecord-marker=4 -ffixed-line-length-0

#PROGRAM = gc2cmaq


#LIBS =    -L/usr/local/ioapi-3.1/Linux2_x86_64 -lioapi -lnetcdf
#INCLUDE = -I/usr/local/ioapi-3.1/ioapi/fixed_src \
          -I/usr/local/netcdf-4.0/include

RAW = \
 $(PROGRAM).f

default:
	@echo '--------------------------------------------------------------------------'
	@echo 'This is the Makefile for gc2cmaq'
	@echo ' '
	@echo 'Syntax is: "make option"'
	@echo ' '
	@echo 'Acceptable options are:'
	@echo '  SAPRC07TIC_AE6I -- for CMAQ SAPRC07TIC gas and AE6I aerosol mechanism from GEOS-Chem'
	@echo '  CB05E51_AE6     -- for CMAQ CB05E51 gas and AE6 aerosol mechanism from GEOS-Chem'
	@echo '--------------------------------------------------------------------------'
	@echo ' '

SAPRC07TIC_AE6I:       
	@rm -f G2Lconv.EXT
	@ln -s G2Lconv_SAPRC07TIC_AE6I_from_GEOS-Chem.EXT G2Lconv.EXT
	$(F90) $(F90_OPTS) $(INCLUDE) -o $(@) $(RAW) $(LIBS)
	@mv SAPRC07TIC_AE6I $(PROGRAM)_SAPRC07TIC_AE6I

CB05E51_AE6:       
	@rm -f G2Lconv.EXT
	@ln -s G2Lconv_CB05E51_AE6_from_GEOS-Chem.EXT G2Lconv.EXT
	$(F90) $(F90_OPTS) $(INCLUDE) -o $(@) $(RAW) $(LIBS)
	@mv CB05E51_AE6 $(PROGRAM)_CB05E51_AE6

$(PROGRAM):	$(RAW)
	$(F90) $(F90_OPTS) $(INCLUDE) -o $(@) $(RAW) $(LIBS)

clean:
	rm -f *.o
