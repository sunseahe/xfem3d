#FC = ifort
FC = gfortran

#CC = icpc
CC = g++

#Folders
src_path = ./src
obj_path = ./obj

#Intel flags fortran
FCFLAGSINTEL = -module $(obj_path) -I$(MKLROOT)/include/intel64/lp64 -I$(MKLROOT)/include -fpp -openmp -standard-semantics
FLFLAGSINTEL = -L$(MKLROOT)/lib/intel64 -L$(IROOT)/lib/intel64 -lmkl_blas95_lp64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread -lm -liomp5
#Intel flags c++
CCFLAGSINTEL = -I$(ABQROOT)/code/include -cxxlib -w -Wno-deprecated -fpermissive -DTYPENAME=typename -D_LINUX_SOURCE -DABQ_LINUX -DABQ_LNX86_64
CLFLAGSINTEL = -cxxlib -Wl,-Bdynamic -L$(ABQROOT)/code/bin -lstandardB -lABQSMAOdbDdbOdb -lABQSMAOdbApi -lABQSMAOdbCore -lABQSMAOdbCoreGeom -lABQSMAOdbAttrEO -lABQSMAAbuBasicUtils -lABQSMABasShared -lABQSMABasCoreUtils -lABQSMAStiCAE_StableTime -lABQSMABasMem -lABQSMAAbuGeom -lABQSMARomDiagEx -lABQSMASspUmaCore -lABQSMASimInterface -lABQSMAMtxCoreModule
#GCC flags fortran
FCFLAGSGCC = -fintrinsic-modules-path $(obj_path) -J$(obj_path) -I$(F95ROOT)/include/intel64/lp64 -I$(MKLROOT)/include -cpp -fopenmp
FLFLAGSGCC = -L$(MKLROOT)/lib/intel64 $(F95ROOT)/lib/intel64/libmkl_blas95_lp64.a $(F95ROOT)/lib/intel64/libmkl_lapack95_lp64.a -lmkl_gf_lp64 -lmkl_core -lmkl_intel_thread -liomp5 -ldl -lpthread -lm
#GCC flags c++
CCFLAGSGCC = -I$(ABQROOT)/code/include -lstdc++ -w -Wno-deprecated -fpermissive -DTYPENAME=typename -D_LINUX_SOURCE -DABQ_LINUX -DABQ_LNX86_64
CLFLAGSGCC =  -lstdc++ -Wl,-Bdynamic -L$(ABQROOT)/code/bin -lstandardB -lABQSMAOdbDdbOdb -lABQSMAOdbApi -lABQSMAOdbCore -lABQSMAOdbCoreGeom -lABQSMAOdbAttrEO -lABQSMAAbuBasicUtils -lABQSMABasShared -lABQSMABasCoreUtils -lABQSMAStiCAE_StableTime -lABQSMABasMem -lABQSMAAbuGeom -lABQSMARomDiagEx -lABQSMASspUmaCore -lABQSMASimInterface -lABQSMAMtxCoreModule

#Source code
vpath %.cpp $(src_path)
vpath %.f90 $(src_path)

SRCF = types.f90 linked_list.f90 determinant.f90 tokenize_string.f90 lsf_test_functions.f90 fe_c3d4.f90 fe_c3d10.f90 xfem_tetra.f90 read_input.f90 write_odb.f90 volume_integral.f90 main.f90
SRCC = odb_routines.cpp

#Objects

OBJCT = $(subst .cpp,.o,$(SRCC))
OBJFT = $(subst .f90,.o,$(SRCF))

OBJC = $(addprefix $(obj_path)/,$(OBJCT))
OBJF = $(addprefix $(obj_path)/,$(OBJFT))

#Intel
ifeq ($(FC),ifort)
  # fortran
  all: FCFLAGS = $(FCFLAGSINTEL)
  all: FLFLAGS = $(FLFLAGSINTEL)
  debug: FCFLAGS += $(FCFLAGSINTEL) -g -debug full -warn all -check all -std08 -diag-error-limit 1 -traceback -DDEBUG
  debug: FLFLAGS = $(FLFLAGSINTEL)
  opt: FCFLAGS += $(FCFLAGSINTEL) -O3
  opt: FLFLAGS = $(FLFLAGSINTEL)
  # c++
  all: CCFLAGS = $(CCFLAGSINTEL)
  all: CLFLAGS = $(CLFLAGSINTEL)
  debug: CCFLAGS += $(CCFLAGSINTEL) -g -Wall -traceback
  debug: CLFLAGS = $(CLFLAGSINTEL)
  opt: CCFLAGS += $(CCFLAGSINTEL) -O3
  opt: CLFLAGS = $(CLFLAGSINTEL)
endif
#GCC
ifeq ($(FC),gfortran)
  # fortran
  all: FCFLAGS = $(FCFLAGSGCC)
  all: FLFLAGS = $(FLFLAGSGCC)
  debug: FCFLAGS += $(FCFLAGSGCC) -g -W -Wall -Wextra -fbounds-check -pedantic -std=f2008 -Wunderflow -O -fbacktrace -ffpe-trap=zero,overflow,underflow -fmax-errors=1 -Wfatal-errors -DDEBUG
  debug: FLFLAGS = $(FLFLAGSGCC)
  opt: FCFLAGS += $(FCFLAGSGCC) -O3 -march=native -ffast-math -funroll-loops
  opt: FLFLAGS = $(FLFLAGSGCC)
  # c++
  all: CCFLAGS = $(CCFLAGSGCC)
  all: CLFLAGS = $(CLFLAGSGCC)
  debug: CCFLAGS += $(CCFLAGSGCC) -g -Wall
  debug: CLFLAGS = $(CLFLAGSGCC)
  opt: CCFLAGS += $(CCFLAGSGCC) -O3
  opt: CLFLAGS = $(CLFLAGSGCC)
endif

PROGRAM  = xfem

all: $(PROGRAM)
debug: $(PROGRAM)
opt: $(PROGRAM)

#COMPILE FORTRAN
$(obj_path)/%.o: %.f90
	@echo $(@F)
	@$(FC) $(FCFLAGS) -c $< -o $(obj_path)/$(@F)
#COMPILE C++
$(obj_path)/%.o: %.cpp
	@echo $(@F)
	@$(CC) $(CCFLAGS) -c $< -o $(obj_path)/$(@F)

#LINK
$(PROGRAM): $(OBJF) $(OBJC)
	@echo 'Linking' $(PROGRAM) '...'
	@$(FC) $^ $(FLFLAGS) $(CLFLAGS) -o $@
	@echo 'Done.'

.PHONY: clean, veryclean

clean:
	rm -f  $(obj_path)/*.o $(obj_path)/*.mod

veryclean:
	rm -f *~ $(PROGRAM)
