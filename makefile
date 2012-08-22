#
# A general makefile with support for several compilers.
#
#-----------------------------------------------------------------------------
# This file is part of BSAM.
#
# BSAM is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# BSAM is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# BSAM.  If not, see <http://www.gnu.org/licenses/>.
#-----------------------------------------------------------------------------
#

# Set some system variables 
#=============================================================================
PROC      = $(shell uname -p)
ifeq ($(PROC),unknown)
  PROC    = $(shell uname -m)
endif
LIB       = ../lib
EXE       = xbsam
OMAIN     = objects
ODIR      = $(OMAIN)/$(TARGET)
BDIR      = ./bin
OUT       = out
BIN       = $(BDIR)/$(EXE)_$(TARGET)
MOD       = -module $(ODIR)
ifeq ($(FOR),gfortran)
  MOD = -J$(ODIR) -I$(ODIR)
else ifeq ($(FOR),sunf90)
  MOD = -moddir=$(ODIR)
else ifeq ($(FOR),g95)
  MOD = -fmod=$(ODIR)
endif

targets   = debug optim
modes     = debug profile optim r16
compilers = ifort gfortran pgf90 g95 sunf90 ifc
default   = gfortran

# Define set of compiler flags for the different compilers
#=============================================================================
# ifort flags
#-----------------------------------------------------------------------------
# Notes:
# * For optim_ifort: -axN ? (is also compatible with generic processor)
flags_ifort_debug   = "-r8 -fpe0 -g -fp-model precise -warn all       \
                       -warn nounused -check all -debug all -traceback"
flags_ifort_profile = "-r8 -fpe0 -g -pg"
flags_ifort_optim   = "-r8 -fpe0 -fast -heap-arrays -w -ax"
flags_ifort_r16     = "-r16 -fpe0 -fast -heap-arrays -w -ax"
#flags_ifort_r16     = "-r16 -fpe0 -g -fp-model precise -warn all -warn unused -traceback -check all -debug all"
ifeq ($(proc),i686)
  flags_ifort_optim = "-r8 -fpe0 -fast -w -xN"
endif

# pgf90 flags
#-----------------------------------------------------------------------------
flags_pgf90_debug   = "-r8 -g -Kieee -Ktrap=inv,divz,ovf -Mbounds \
		      -Mchkfpstk -Mchkptr -Minform=inform"
flags_pgf90_profile = "-r8 -pg -O0 -Mprof=lines"
flags_pgf90_optim   = "-r8 -fastsse -Mipa=fast,safe -Ktrap=inv,divz,ovf"

# gfortran flags
#-----------------------------------------------------------------------------
# Optional flags, may be useful for other architectures
#  -mtune=nocona
#  -cpu=970
#  -mtune=970
#  -powerpc64
#  -mpowerpc-gpopt
#  -force-cpusubtype_ALL
#  -ffast-math 
flags_gfortran_debug   = "-fdefault-real-8 -g -fbounds-check \
                          -ffpe-trap=invalid,zero,overflow -Wall -Wextra \
                          -Wno-unused-parameter -Wno-unused-dummy-argument"
flags_gfortran_profile = "-fdefault-real-8 -g -pg"
flags_gfortran_optim   = "-fdefault-real-8 -O3 -funroll-loops \
			  -march=native -msse3"

# sunf90 flags
#-----------------------------------------------------------------------------
flags_sunf90_debug  = "-xtypemap=real:64 -C -g -w4 -stackvar"
flags_sunf90_optim  = "-xtypemap=real:64 -fast"

# g95 flags
#-----------------------------------------------------------------------------
# Other options that might be useful
#  -ffast-math
#  -march=nocona
#
flags_g95_debug  = "-r8 -g -fbounds-check -ftrace=full -mieee-fp -Wall \
                    -fno-second-underscore -Wextra"
flags_g95_optim  = "-r8 -O3 -msse3 -funroll-loops -fno-second-underscore"
ifneq ($(proc),x86_64)
  flags_g95_optim += "-malign-double"
endif

# ifc flags
#-----------------------------------------------------------------------------
flags_ifc_optim = "ifc -i_dynamic -O3 -tpp7 -I$(ODIR) -module $(ODIR)"
flags_ifc_optim = "ifc -i_dynamic -CB -I$(ODIR) -module $(ODIR)"

# Define some special targets
#=============================================================================
.PHONY: help clean dist_clean build xbsam

xbsam: optim
	./xbsam

help:
	@echo "Usage: make <target>"
	@echo "<target> is one of the following:"
	@echo
	@echo "  clean"
	@echo "  help"
	@echo " $(foreach t,$(sort $(targets)), $(t)\n)"
	@echo "The default compiler target is $(default)."

clean:
	rm -f $(EXE)
	rm -f $(BDIR)/*
	rm -f $(OMAIN)/*/*

dist_clean: clean
	rm -rf $(OMAIN)
	rm -rf $(BDIR)
	rm -rf $(OUT)

build:
	rm -f $(EXE)
	mkdir -p $(ODIR)
	mkdir -p $(BDIR)
	mkdir -p $(OUT)
	$(MAKE) -e $(BIN)
	ln -s $(BIN) $(EXE)

# Set default targets
#=============================================================================
debug: debug_$(default)

optim: optim_$(default)

# Create the compiler specific targets
#=============================================================================

# Define objects
OBJ = $(ODIR)/nodeinfodef.o     \
      $(ODIR)/problemdef.o      \
      $(ODIR)/treeops.o         \
      $(ODIR)/gridutilities.o   \
      $(ODIR)/problem.o         \
      $(ODIR)/bsamstorage.o     \
      $(ODIR)/boundary.o        \
      $(ODIR)/bsaminputoutput.o \
      $(ODIR)/afasroutines.o    \
      $(ODIR)/bsamroutines.o    \
      $(ODIR)/bsamdriver.o

# Use a macro (or canned sequence) to automatically create the targets
define create_target
  ifdef flags_$(2)_$(1)
    targets += $(1)_$(2)
$(1)_$(2):
	$(MAKE) FOR=$(2) FFLAGS=$(flags_$(2)_$(1)) TARGET=$$@ -e build
  endif
endef
$(foreach mode,$(modes),$(foreach comp,$(compilers),\
  $(eval $(call create_target,$(mode),$(comp)))))

# Create targets for the objects
make_macro = $(FOR) $(FFLAGS) $(MOD) -c $(1) -o $(2)
$(ODIR)/nodeinfodef.o: $(LIB)/nodeinfodef.f90
	$(FOR) $(FFLAGS) $(MOD) -o $@ -c $(LIB)/nodeinfodef.f90
$(ODIR)/problemdef.o: problemdef.f90
	$(FOR) $(FFLAGS) $(MOD) -o $@ -c problemdef.f90
$(ODIR)/treeops.o: $(LIB)/treeops.f90
	$(FOR) $(FFLAGS) $(MOD) -o $@ -c $(LIB)/treeops.f90
$(ODIR)/gridutilities.o: $(LIB)/gridutilities.f90
	$(FOR) $(FFLAGS) $(MOD) -o $@ -c $(LIB)/gridutilities.f90
$(ODIR)/problem.o: problem.f90
	$(FOR) $(FFLAGS) $(MOD) -o $@ -c problem.f90
$(ODIR)/bsamstorage.o: $(LIB)/bsamstorage.f90
	$(FOR) $(FFLAGS) $(MOD) -o $@ -c $(LIB)/bsamstorage.f90
$(ODIR)/boundary.o: $(LIB)/boundary.f90
	$(FOR) $(FFLAGS) $(MOD) -o $@ -c $(LIB)/boundary.f90
$(ODIR)/bsaminputoutput.o: $(LIB)/bsaminputoutput.f90
	$(FOR) $(FFLAGS) $(MOD) -o $@ -c $(LIB)/bsaminputoutput.f90
$(ODIR)/afasroutines.o: $(LIB)/afasroutines.f90
	$(FOR) $(FFLAGS) $(MOD) -o $@ -c $(LIB)/afasroutines.f90
$(ODIR)/bsamroutines.o: $(LIB)/bsamroutines.f90
	$(FOR) $(FFLAGS) $(MOD) -o $@ -c $(LIB)/bsamroutines.f90
$(ODIR)/bsamdriver.o: $(LIB)/bsamdriver.f90
	$(FOR) $(FFLAGS) $(MOD) -o $@ -c $(LIB)/bsamdriver.f90

# Create the target for the executable
$(BIN): $(OBJ)
	$(FOR) $(FFLAGS) -o $(BIN) $(OBJ)
