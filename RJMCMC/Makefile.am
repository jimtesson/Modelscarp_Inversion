
ACLOCAL_AMFLAGS = -I m4

SUBDIRS = src \
	include

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA= rjmcmc.pc

DIST_SUBDIRS = src \
	include \
	tests \
	examples/1d/partitioned/fm/functionfit \
	examples/1d/partitioned/fm/functionfitf \
	examples/1d/partitioned/fm/gaussian \
	examples/1d/partitioned/fm/natural \
	examples/1d/partitioned/fm/regression \
	examples/1d/partitioned/fm/step \
	examples/1d/partitioned/fm/synthdata \
	examples/1d/partitioned/fm/velocityf \
	examples/1d/partitioned/regression/multiquad \
	examples/1d/partitioned/regression/multistep \
	examples/1d/partitioned/regression/sawtooth \
	examples/1d/partitioned/regression/zeromultistep \
	examples/1d/single/regression/cubic \
	examples/1d/single/fm/gaussian \
	examples/1d/single/fm/simplef \
	examples/1d/single/fm/simpleimage \
	examples/1d/single/fm/spherefit \
	examples/2d/partitioned/fm/mpi_regression \
	examples/2d/partitioned/fm/regression \
	examples/2d/partitioned/fm/gaussian \
	examples/2d/partitioned/regression/ausmoho \
	examples/2d/partitioned/regression/callbackex \
	examples/2d/partitioned/regression/disc \
	examples/2d/partitioned/regression/gaussian \
	examples/2d/partitioned/regression/square \
	tests/regression/1d/partitioned/fm/MPI_part1d_forwardmodel_c \
	tests/regression/1d/partitioned/fm/MPI_part1d_forwardmodel_hierarchical_c \
	tests/regression/1d/partitioned/fm/MPI_part1d_forwardmodel_natural_c \
	tests/regression/1d/partitioned/fm/MPI_part1d_forwardmodel_natural_hierarchical_c \
	tests/regression/1d/partitioned/fm/part1d_forwardmodel_c \
	tests/regression/1d/partitioned/fm/part1d_forwardmodel_hierarchical_c \
	tests/regression/1d/partitioned/fm/part1d_forwardmodel_natural_c \
	tests/regression/1d/partitioned/fm/part1d_forwardmodel_natural_hierarchical_c \
	tests/regression/1d/partitioned/fm/part1d_forwardmodel_zero_cubic_c \
	tests/regression/1d/single/fm/MPI_single_forwardmodel_c \
	tests/regression/1d/single/fm/MPI_single_forwardmodel_hierarchical_c \
	tests/regression/1d/single/fm/single_forwardmodel_c \
	tests/regression/1d/single/fm/single_forwardmodel_hierarchical_c \
	tests/regression/2d/partitioned/fm/MPI_part2d_forwardmodel_c \
	tests/regression/2d/partitioned/fm/MPI_part2d_forwardmodel_hierarchical_c \
	tests/regression/2d/partitioned/fm/part2d_forwardmodel_c \
	tests/regression/2d/partitioned/fm/part2d_forwardmodel_hierarchical_c

EXTRA_DIST = \
	macosx/builduniversal.sh \
	macosx/installuniversal.sh \
	macosx/buildforpython.sh \
	doc/rjmcmc.dox \
	doc/images/main_tree.py \
	doc/images/Makefile \
	python/swig/rjmcmc.i \
	python/swig/rjmcmc_helper.c \
	python/swig/rjmcmc_helper.h \
	python/swig/setup.py \
	python/tutorial/fm/ch0-createdata.py \
	python/tutorial/fm/ch1-loading.py \
	python/tutorial/fm/ch2-analyse.py \
	python/tutorial/fm/ch3-regression.py \
	python/tutorial/single/data.txt \
	python/tutorial/single/ch0-createdata.py \
	python/tutorial/single/ch1-loading.py \
	python/tutorial/single/ch2-analyse.py \
	python/tutorial/single/ch3-orderanalysis.py \
	python/tutorial/single/ch4-confidence.py \
	python/tutorial/single/ch5-hierarchical.py \
	python/tutorial/single/doc/tufte-book.cls \
	python/tutorial/single/doc/tufte-common.def \
	python/tutorial/single/doc/tufte-handout.cls \
	python/tutorial/single/doc/tutorial_single.tex \
	python/tutorial/multi/data.txt \
	python/tutorial/multi/ch0-createdata.py \
	python/tutorial/multi/ch1-loading.py \
	python/tutorial/multi/ch2-analyse.py \
	python/tutorial/multi/ch3-orderanalysis.py \
	python/tutorial/multi/ch4-confidence.py \
	python/tutorial/multi/ch5-hierarchical.py \
	python/tutorial/multi/doc/tufte-book.cls \
	python/tutorial/multi/doc/tufte-common.def \
	python/tutorial/multi/doc/tufte-handout.cls \
	python/tutorial/multi/doc/tutorial_multi.tex \
	doc/images/main_tree.pdf \
	doc/images/main_tree.png \
	template/part1d/Makefile \
	template/part1d/part1dtemplate.c \
	template/part1d/data.txt \
	template/part1d_fortran/Makefile \
	template/part1d_fortran/part1d_fortran.f90 \
	template/part1d_fortran/data.txt \
	template/part1d_hierarchical_fortran/create_data.py \
	template/part1d_hierarchical_fortran/data.txt \
	template/part1d_hierarchical_fortran/Makefile \
	template/part1d_hierarchical_fortran/part1d_hierarchical_fortran.f90 \
	template/part1d_hierarchical_fortran/realdata.txt \
	template/part1d_hierarchical_fortran/svdbacksub.f90 \
	template/part1d_hierarchical_fortran/svdtest.f90 \
	template/part1d_hierarchical_fortran/svdtest.m \
	template/part1d_hierarchical_fortran/blas/caxpy.f     \
	template/part1d_hierarchical_fortran/blas/ccopy.f	 \
	template/part1d_hierarchical_fortran/blas/cdotc.f	 \
	template/part1d_hierarchical_fortran/blas/cdotu.f	 \
	template/part1d_hierarchical_fortran/blas/cgbmv.f	 \
	template/part1d_hierarchical_fortran/blas/cgemm.f	 \
	template/part1d_hierarchical_fortran/blas/cgemv.f	 \
	template/part1d_hierarchical_fortran/blas/cgerc.f	 \
	template/part1d_hierarchical_fortran/blas/cgeru.f	 \
	template/part1d_hierarchical_fortran/blas/chbmv.f	 \
	template/part1d_hierarchical_fortran/blas/chemm.f	 \
	template/part1d_hierarchical_fortran/blas/chemv.f	 \
	template/part1d_hierarchical_fortran/blas/cher2.f	 \
	template/part1d_hierarchical_fortran/blas/cher2k.f	 \
	template/part1d_hierarchical_fortran/blas/cher.f	 \
	template/part1d_hierarchical_fortran/blas/cherk.f	 \
	template/part1d_hierarchical_fortran/blas/chpmv.f	 \
	template/part1d_hierarchical_fortran/blas/chpr2.f	 \
	template/part1d_hierarchical_fortran/blas/chpr.f	 \
	template/part1d_hierarchical_fortran/blas/crotg.f	 \
	template/part1d_hierarchical_fortran/blas/cscal.f	 \
	template/part1d_hierarchical_fortran/blas/csrot.f	 \
	template/part1d_hierarchical_fortran/blas/csscal.f	 \
	template/part1d_hierarchical_fortran/blas/cswap.f	 \
	template/part1d_hierarchical_fortran/blas/csymm.f	 \
	template/part1d_hierarchical_fortran/blas/csyr2k.f	 \
	template/part1d_hierarchical_fortran/blas/csyrk.f	 \
	template/part1d_hierarchical_fortran/blas/ctbmv.f	 \
	template/part1d_hierarchical_fortran/blas/ctbsv.f	 \
	template/part1d_hierarchical_fortran/blas/ctpmv.f	 \
	template/part1d_hierarchical_fortran/blas/ctpsv.f	 \
	template/part1d_hierarchical_fortran/blas/ctrmm.f	 \
	template/part1d_hierarchical_fortran/blas/ctrmv.f	 \
	template/part1d_hierarchical_fortran/blas/ctrsm.f	 \
	template/part1d_hierarchical_fortran/blas/ctrsv.f	 \
	template/part1d_hierarchical_fortran/blas/dasum.f	 \
	template/part1d_hierarchical_fortran/blas/daxpy.f	 \
	template/part1d_hierarchical_fortran/blas/dcabs1.f	 \
	template/part1d_hierarchical_fortran/blas/dcopy.f	 \
	template/part1d_hierarchical_fortran/blas/ddot.f	 \
	template/part1d_hierarchical_fortran/blas/dgbmv.f	 \
	template/part1d_hierarchical_fortran/blas/dgemm.f	 \
	template/part1d_hierarchical_fortran/blas/dgemv.f	 \
	template/part1d_hierarchical_fortran/blas/dger.f	 \
	template/part1d_hierarchical_fortran/blas/dnrm2.f	 \
	template/part1d_hierarchical_fortran/blas/drot.f	 \
	template/part1d_hierarchical_fortran/blas/drotg.f	 \
	template/part1d_hierarchical_fortran/blas/drotm.f	 \
	template/part1d_hierarchical_fortran/blas/drotmg.f	 \
	template/part1d_hierarchical_fortran/blas/dsbmv.f	 \
	template/part1d_hierarchical_fortran/blas/dscal.f	 \
	template/part1d_hierarchical_fortran/blas/dsdot.f	 \
	template/part1d_hierarchical_fortran/blas/dspmv.f	 \
	template/part1d_hierarchical_fortran/blas/dspr2.f	 \
	template/part1d_hierarchical_fortran/blas/dspr.f	 \
	template/part1d_hierarchical_fortran/blas/dswap.f	 \
	template/part1d_hierarchical_fortran/blas/dsymm.f	 \
	template/part1d_hierarchical_fortran/blas/dsymv.f	 \
	template/part1d_hierarchical_fortran/blas/dsyr2.f	 \
	template/part1d_hierarchical_fortran/blas/dsyr2k.f	 \
	template/part1d_hierarchical_fortran/blas/dsyr.f	 \
	template/part1d_hierarchical_fortran/blas/dsyrk.f	 \
	template/part1d_hierarchical_fortran/blas/dtbmv.f	 \
	template/part1d_hierarchical_fortran/blas/dtbsv.f	 \
	template/part1d_hierarchical_fortran/blas/dtpmv.f	 \
	template/part1d_hierarchical_fortran/blas/dtpsv.f	 \
	template/part1d_hierarchical_fortran/blas/dtrmm.f	 \
	template/part1d_hierarchical_fortran/blas/dtrmv.f	 \
	template/part1d_hierarchical_fortran/blas/dtrsm.f	 \
	template/part1d_hierarchical_fortran/blas/dtrsv.f	 \
	template/part1d_hierarchical_fortran/blas/dzasum.f	 \
	template/part1d_hierarchical_fortran/blas/dznrm2.f	 \
	template/part1d_hierarchical_fortran/blas/icamax.f	 \
	template/part1d_hierarchical_fortran/blas/idamax.f	 \
	template/part1d_hierarchical_fortran/blas/isamax.f	 \
	template/part1d_hierarchical_fortran/blas/izamax.f	 \
	template/part1d_hierarchical_fortran/blas/lsame.f	 \
	template/part1d_hierarchical_fortran/blas/Makefile	 \
	template/part1d_hierarchical_fortran/blas/sasum.f	 \
	template/part1d_hierarchical_fortran/blas/saxpy.f	 \
	template/part1d_hierarchical_fortran/blas/scabs1.f	 \
	template/part1d_hierarchical_fortran/blas/scasum.f	 \
	template/part1d_hierarchical_fortran/blas/scnrm2.f	 \
	template/part1d_hierarchical_fortran/blas/scopy.f	 \
	template/part1d_hierarchical_fortran/blas/sdot.f	 \
	template/part1d_hierarchical_fortran/blas/sdsdot.f	 \
	template/part1d_hierarchical_fortran/blas/sgbmv.f	 \
	template/part1d_hierarchical_fortran/blas/sgemm.f	 \
	template/part1d_hierarchical_fortran/blas/sgemv.f	 \
	template/part1d_hierarchical_fortran/blas/sger.f	 \
	template/part1d_hierarchical_fortran/blas/snrm2.f	 \
	template/part1d_hierarchical_fortran/blas/srot.f	 \
	template/part1d_hierarchical_fortran/blas/srotg.f	 \
	template/part1d_hierarchical_fortran/blas/srotm.f	 \
	template/part1d_hierarchical_fortran/blas/srotmg.f	 \
	template/part1d_hierarchical_fortran/blas/ssbmv.f	 \
	template/part1d_hierarchical_fortran/blas/sscal.f	 \
	template/part1d_hierarchical_fortran/blas/sspmv.f	 \
	template/part1d_hierarchical_fortran/blas/sspr2.f	 \
	template/part1d_hierarchical_fortran/blas/sspr.f	 \
	template/part1d_hierarchical_fortran/blas/sswap.f	 \
	template/part1d_hierarchical_fortran/blas/ssymm.f	 \
	template/part1d_hierarchical_fortran/blas/ssymv.f	 \
	template/part1d_hierarchical_fortran/blas/ssyr2.f	 \
	template/part1d_hierarchical_fortran/blas/ssyr2k.f	 \
	template/part1d_hierarchical_fortran/blas/ssyr.f	 \
	template/part1d_hierarchical_fortran/blas/ssyrk.f	 \
	template/part1d_hierarchical_fortran/blas/stbmv.f	 \
	template/part1d_hierarchical_fortran/blas/stbsv.f	 \
	template/part1d_hierarchical_fortran/blas/stpmv.f	 \
	template/part1d_hierarchical_fortran/blas/stpsv.f	 \
	template/part1d_hierarchical_fortran/blas/strmm.f	 \
	template/part1d_hierarchical_fortran/blas/strmv.f	 \
	template/part1d_hierarchical_fortran/blas/strsm.f	 \
	template/part1d_hierarchical_fortran/blas/strsv.f	 \
	template/part1d_hierarchical_fortran/blas/xerbla_array.f \
	template/part1d_hierarchical_fortran/blas/xerbla.f	 \
	template/part1d_hierarchical_fortran/blas/zaxpy.f	 \
	template/part1d_hierarchical_fortran/blas/zcopy.f	 \
	template/part1d_hierarchical_fortran/blas/zdotc.f	 \
	template/part1d_hierarchical_fortran/blas/zdotu.f	 \
	template/part1d_hierarchical_fortran/blas/zdrot.f	 \
	template/part1d_hierarchical_fortran/blas/zdscal.f	 \
	template/part1d_hierarchical_fortran/blas/zgbmv.f	 \
	template/part1d_hierarchical_fortran/blas/zgemm.f	 \
	template/part1d_hierarchical_fortran/blas/zgemv.f	 \
	template/part1d_hierarchical_fortran/blas/zgerc.f	 \
	template/part1d_hierarchical_fortran/blas/zgeru.f	 \
	template/part1d_hierarchical_fortran/blas/zhbmv.f	 \
	template/part1d_hierarchical_fortran/blas/zhemm.f	 \
	template/part1d_hierarchical_fortran/blas/zhemv.f	 \
	template/part1d_hierarchical_fortran/blas/zher2.f	 \
	template/part1d_hierarchical_fortran/blas/zher2k.f	 \
	template/part1d_hierarchical_fortran/blas/zher.f	 \
	template/part1d_hierarchical_fortran/blas/zherk.f	 \
	template/part1d_hierarchical_fortran/blas/zhpmv.f	 \
	template/part1d_hierarchical_fortran/blas/zhpr2.f	 \
	template/part1d_hierarchical_fortran/blas/zhpr.f	 \
	template/part1d_hierarchical_fortran/blas/zrotg.f	 \
	template/part1d_hierarchical_fortran/blas/zscal.f	 \
	template/part1d_hierarchical_fortran/blas/zswap.f	 \
	template/part1d_hierarchical_fortran/blas/zsymm.f	 \
	template/part1d_hierarchical_fortran/blas/zsyr2k.f	 \
	template/part1d_hierarchical_fortran/blas/zsyrk.f	 \
	template/part1d_hierarchical_fortran/blas/ztbmv.f	 \
	template/part1d_hierarchical_fortran/blas/ztbsv.f	 \
	template/part1d_hierarchical_fortran/blas/ztpmv.f	 \
	template/part1d_hierarchical_fortran/blas/ztpsv.f	 \
	template/part1d_hierarchical_fortran/blas/ztrmm.f	 \
	template/part1d_hierarchical_fortran/blas/ztrmv.f	 \
	template/part1d_hierarchical_fortran/blas/ztrsm.f	 \
	template/part1d_hierarchical_fortran/blas/ztrsv.f     \
	template/part1d_hierarchical_fortran/lapack/dbdsdc.f   \
	template/part1d_hierarchical_fortran/lapack/dbdsqr.f   \
	template/part1d_hierarchical_fortran/lapack/dgebd2.f   \
	template/part1d_hierarchical_fortran/lapack/dgebrd.f   \
	template/part1d_hierarchical_fortran/lapack/dgelq2.f   \
	template/part1d_hierarchical_fortran/lapack/dgelqf.f   \
	template/part1d_hierarchical_fortran/lapack/dgeqr2.f   \
	template/part1d_hierarchical_fortran/lapack/dgeqrf.f   \
	template/part1d_hierarchical_fortran/lapack/dgesdd.f   \
	template/part1d_hierarchical_fortran/lapack/disnan.f   \
	template/part1d_hierarchical_fortran/lapack/dlabrd.f   \
	template/part1d_hierarchical_fortran/lapack/dlacpy.f   \
	template/part1d_hierarchical_fortran/lapack/dlaed6.f   \
	template/part1d_hierarchical_fortran/lapack/dlaisnan.f \
	template/part1d_hierarchical_fortran/lapack/dlamch.f   \
	template/part1d_hierarchical_fortran/lapack/dlamrg.f   \
	template/part1d_hierarchical_fortran/lapack/dlange.f   \
	template/part1d_hierarchical_fortran/lapack/dlanst.f   \
	template/part1d_hierarchical_fortran/lapack/dlapy2.f   \
	template/part1d_hierarchical_fortran/lapack/dlarfb.f   \
	template/part1d_hierarchical_fortran/lapack/dlarf.f    \
	template/part1d_hierarchical_fortran/lapack/dlarfg.f   \
	template/part1d_hierarchical_fortran/lapack/dlarft.f   \
	template/part1d_hierarchical_fortran/lapack/dlartg.f   \
	template/part1d_hierarchical_fortran/lapack/dlas2.f    \
	template/part1d_hierarchical_fortran/lapack/dlascl.f   \
	template/part1d_hierarchical_fortran/lapack/dlasd0.f   \
	template/part1d_hierarchical_fortran/lapack/dlasd1.f   \
	template/part1d_hierarchical_fortran/lapack/dlasd2.f   \
	template/part1d_hierarchical_fortran/lapack/dlasd3.f   \
	template/part1d_hierarchical_fortran/lapack/dlasd4.f   \
	template/part1d_hierarchical_fortran/lapack/dlasd5.f   \
	template/part1d_hierarchical_fortran/lapack/dlasd6.f   \
	template/part1d_hierarchical_fortran/lapack/dlasd7.f   \
	template/part1d_hierarchical_fortran/lapack/dlasd8.f   \
	template/part1d_hierarchical_fortran/lapack/dlasda.f   \
	template/part1d_hierarchical_fortran/lapack/dlasdq.f   \
	template/part1d_hierarchical_fortran/lapack/dlasdt.f   \
	template/part1d_hierarchical_fortran/lapack/dlaset.f   \
	template/part1d_hierarchical_fortran/lapack/dlasq1.f   \
	template/part1d_hierarchical_fortran/lapack/dlasq2.f   \
	template/part1d_hierarchical_fortran/lapack/dlasq3.f   \
	template/part1d_hierarchical_fortran/lapack/dlasq4.f   \
	template/part1d_hierarchical_fortran/lapack/dlasq5.f   \
	template/part1d_hierarchical_fortran/lapack/dlasq6.f   \
	template/part1d_hierarchical_fortran/lapack/dlasr.f    \
	template/part1d_hierarchical_fortran/lapack/dlasrt.f   \
	template/part1d_hierarchical_fortran/lapack/dlassq.f   \
	template/part1d_hierarchical_fortran/lapack/dlasv2.f   \
	template/part1d_hierarchical_fortran/lapack/dorg2r.f   \
	template/part1d_hierarchical_fortran/lapack/dorgbr.f   \
	template/part1d_hierarchical_fortran/lapack/dorgl2.f   \
	template/part1d_hierarchical_fortran/lapack/dorglq.f   \
	template/part1d_hierarchical_fortran/lapack/dorgqr.f   \
	template/part1d_hierarchical_fortran/lapack/dorm2r.f   \
	template/part1d_hierarchical_fortran/lapack/dormbr.f   \
	template/part1d_hierarchical_fortran/lapack/dorml2.f   \
	template/part1d_hierarchical_fortran/lapack/dormlq.f   \
	template/part1d_hierarchical_fortran/lapack/dormqr.f   \
	template/part1d_hierarchical_fortran/lapack/ieeeck.f   \
	template/part1d_hierarchical_fortran/lapack/iladlc.f   \
	template/part1d_hierarchical_fortran/lapack/iladlr.f   \
	template/part1d_hierarchical_fortran/lapack/ilaenv.f   \
	template/part1d_hierarchical_fortran/lapack/iparmq.f   \
	template/part1d_hierarchical_fortran/lapack/lsame.f    \
	template/part1d_hierarchical_fortran/lapack/Makefile   \
	template/part1d_hierarchical_fortran/lapack/xerbla.f 









