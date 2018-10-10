
include $(top_srcdir)/Makefile.am.common

SUBDIRS = .

# fastasplit can clash with fastaplit from exonerate, so we call it fastasplit_random
bin_PROGRAMS = countsnps fastasplit_random indicate translate

PROG_LIBS = \
$(NULL)

countsnps_SOURCES = countsnps.c
fastasplit_random_SOURCES = fastasplit.c
indicate_SOURCES = indicate.c indicate.h
translate_SOURCES = translate_main.c

countsnps_LDADD = $(PROG_LIBS)
fastasplit_random_LDADD = $(PROG_LIBS)
translate_LDADD = $(PROG_LIBS)

indicate_CPPFLAGS = -I../OBDA
indicate_LDADD = $(PROG_LIBS)
indicate_LDADD += $(top_builddir)/OBDA/libOBDA.la


if HAVE_MYSQL_LIB
if HAVE_SAMTOOLS
  bin_PROGRAMS += bamcount bamcount_exon bamcov bammap

  PROG_LIBS += -lhts -lz

  bamcount_SOURCES = bamcount.c bamhelper.h
  bamcount_exon_SOURCES = bamcount_exon.c bamhelper.h
  bamcov_SOURCES = bamcov.c
  bammap_SOURCES = bammap.c bamhelper.h

  bamcount_LDADD = $(PROG_LIBS)
  bamcount_exon_LDADD = $(PROG_LIBS)
  bamcov_LDADD = $(PROG_LIBS)
  bammap_LDADD = $(PROG_LIBS)

if HAVE_LIBCONFIG
  bin_PROGRAMS += RefineSolexaGenes
  RefineSolexaGenes_SOURCES = RefineSolexaGenes.c RefineSolexaGenes.h 

  PROG_LIBS += -lconfig $(MYSQL_LDFLAGS)
  RefineSolexaGenes_LDADD = $(PROG_LIBS)
if HAVE_LIBTCMALLOC
  PROG_LIBS += -ltcmalloc
endif
endif
endif
endif



# Extra files to remove for the maintainer-clean target.
#
MAINTAINERCLEANFILES = $(top_srcdir)/Programs/Makefile.in

