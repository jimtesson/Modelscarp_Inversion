
bin_PROGRAMS = gaussian

gaussian_SOURCES = gaussian.c
gaussian_CPPFLAGS = -I../../../../../include
gaussian_LDADD = ../../../../../src/.libs/librjmcmc.a

if WITHOPENMPI
gaussian_LDADD += -lmpi
endif

EXTRA_DIST = README \
	mkdata.py \
	data.txt \
	plot_mean.py \
	plot_median.py \
	plot_partition_count_histogram.py \
	plot.py


