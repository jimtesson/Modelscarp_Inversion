
ACLOCAL_AMFLAGS = -I m4


MODEL_SOURCES = src/modelscarp/aicc.f \
	src/modelscarp/clcoll.f \
	src/modelscarp/clrock.f src/modelscarp/d_function_2.f \
	src/modelscarp/d_function.f src/modelscarp/fitexp.f \
	src/modelscarp/nelmin.f src/modelscarp/fn.f \
	src/modelscarp/forward.f src/modelscarp/scdepth.f \
	src/modelscarp/scrock.f src/modelscarp/scsurf.f \
	src/modelscarp/sd_function.f src/modelscarp/Ss_function.f

SHARED_SOURCES = src/shared/rf_types.f90 \
	src/shared/utilities.f90 \
	src/shared/forward_model.f90

bin_PROGRAMS =$(top_builddir)/bin/Modelscarp_inv

__top_builddir__bin_Modelscarp_inv_SOURCES = $(MODEL_SOURCES) $(SHARED_SOURCES)  src/Modelscarp_Inv_mpi.f90 
__top_builddir__bin_Modelscarp_inv_FCFLAGS = $(RJMCMC_CFLAGS) -I/usr/include
__top_builddir__bin_Modelscarp_inv_LDADD = $(RJMCMC_LIBS) 
__top_builddir__bin_Modelscarp_inv_DEPENDENCIES =




