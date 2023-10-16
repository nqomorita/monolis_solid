#> Makefile

##> compiler setting
FC     = mpif90
FFLAGS = -fPIC -O2 -mtune=native -march=native -std=legacy -Wno-missing-include-dirs
CC     = mpicc -std=c99
CFLAGS = -fPIC -O2

##> directory setting
# metis library
METIS_DIR  = ./submodule/monolis
METIS_INC  = -I $(METIS_DIR)/include
METIS_LIB  = -L$(METIS_DIR)/lib -lmetis

# mumps library
#MUMPS_DIR  = ./
#MUMPS_INC  = -I $(MUMPS_DIR)/include
#MUMPS_LIB  = -L$(MUMPS_DIR)/lib -lpord -lmumps_common -ldmumps
#-L/usr/local/lib -lscalapack -L/usr/local/Cellar/openblas/0.3.13/lib  -lopenblas

# monolis library
MONOLIS_DIR = ./submodule/monolis
MONOLIS_INC = -I $(MONOLIS_DIR)/include
MONOLIS_LIB = -L$(MONOLIS_DIR)/lib -lmonolis_solver -lgedatsu -lmonolis_utils

LIBS     = $(MONOLIS_LIB) $(MUMPS_LIB) $(METIS_LIB)
INCLUDE  = -I ./include $(MONOLIS_INC)
MOD_DIR  = -J ./include
BIN_DIR  = ./bin
SRC_DIR  = ./src
OBJ_DIR  = ./obj

##> other commands
MAKE = make
CD   = cd
CP   = cp
RM   = rm -rf
AR   = - ar ruv

##> **********
##> target
BIN1     = monolis_solid_l_static
BIN2     = monolis_solid_nl_static
BIN3     = monolis_solid_l_dynamic
BIN4     = monolis_solid_nl_dynamic
TARGET1  = $(addprefix $(BIN_DIR)/, $(BIN1))
TARGET2  = $(addprefix $(BIN_DIR)/, $(BIN2))
TARGET3  = $(addprefix $(BIN_DIR)/, $(BIN3))
TARGET4  = $(addprefix $(BIN_DIR)/, $(BIN4))

SRC_LIST = \
sys/util.f90 \
io/debug.f90 \
io/io.f90 \
sys/solver.f90 \
material/el.f90 \
material/elpl.f90 \
elem/element_C3D8.f90 \
solid/matrix.f90 \
solid/update.f90

SOURCES1  = $(addprefix $(SRC_DIR)/, $(SRC_LIST)) ./src/solid_l_static.f90
SOURCES2  = $(addprefix $(SRC_DIR)/, $(SRC_LIST)) ./src/solid_nl_static.f90
SOURCES3  = $(addprefix $(SRC_DIR)/, $(SRC_LIST)) ./src/solid_l_dynamic.f90
SOURCES4  = $(addprefix $(SRC_DIR)/, $(SRC_LIST)) ./src/solid_nl_dynamic.f90

OBJS1     = $(subst $(SRC_DIR), $(OBJ_DIR), $(SOURCES1:.f90=.o))
OBJS2     = $(subst $(SRC_DIR), $(OBJ_DIR), $(SOURCES2:.f90=.o))
OBJS3     = $(subst $(SRC_DIR), $(OBJ_DIR), $(SOURCES3:.f90=.o))
OBJS4     = $(subst $(SRC_DIR), $(OBJ_DIR), $(SOURCES4:.f90=.o))

all: \
	$(TARGET1)

$(TARGET1): $(OBJS1)
	$(FC) -o $@ $(OBJS1) $(LIBS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FFLAGS) $(INCLUDE) $(MOD_DIR) -o $@ -c $<

clean:
	$(RM) \
	$(OBJS1) \
	$(OBJS2) \
	$(OBJS3) \
	$(OBJS4) \
	$(TARGET1) \
	$(TARGET2) \
	$(TARGET3) \
	$(TARGET4) \
	./include/*.mod

.PHONY: clean
