#
# The main EnsC makefile.am
# Process this file with automake to create Makefile.in.
#

SUBDIRS = OBDA Programs

ACLOCAL_AMFLAGS = -I m4

# keep libtool up to date automatically.
LIBTOOL_DEPS = @LIBTOOL_DEPS@
libtool: $(LIBTOOL_DEPS)
	$(SHELL) ./config.status libtool


# Extra files to remove for the maintainer-clean target.
# Note you cannot use this target to remove directories,
# hence the extra "local" target.
#
MAINTAINERCLEANFILES = \
$(top_srcdir)/Makefile.in \
$(top_srcdir)/aclocal.m4 \
$(top_srcdir)/config.h.in \
$(top_srcdir)/config.h.in~ \
$(top_srcdir)/configure \
$(top_srcdir)/ensc_annotation_tools_version.m4

GENERATED_FILES = $(top_srcdir)/config $(top_srcdir)/autom4te.cache

maintainer-clean-local:
	rm -rf $(GENERATED_FILES)
