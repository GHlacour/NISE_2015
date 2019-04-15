# Makefile
MACHINE=peri

### My settings
ifeq ($(MACHINE), mymachine)
CC=        # Your compiler name
LAPACK=    # Libraries for lapack
LAPACKL=   # Location of liblapack
FFTW3=     # Libraries for fftw3
FFTW3L=    # Location of libfftw3
FFTW3I=    # Location of fftw3.h file
CC_FLAGS=  # Optional compiler flags
endif

### My settings
ifeq ($(MACHINE), peri)
CC = icc
FFTW3=-lfftw3
FFTW3L = -L/cm/shared/apps/intel/compilers_and_libraries/linux/lib/intel64
FFTW3I = -qopenmp
CC_FLAGS = -B -O3 -lm -qopenmp
LAPACK = -llapack -lblas -lpthread
endif

### My settings
ifeq ($(MACHINE), osx9)
CC=gcc
LAPACK=-framework Accelerate
LAPACKL=
FFTW3= -lfftw3
FFTW3L=-L/usr/local/lib
FFTW3I=-I/usr/local/include
CC_FLAGS=-O3 -lm
endif

### My settings
ifeq ($(MACHINE), osx)
CC=gcc
LAPACK=-framework veclib    
LAPACKL=
FFTW3= -lfftw3
FFTW3L= -L/usr/local/lib   
FFTW3I=    
CC_FLAGS= -O3 -lm
endif

### My settings
ifeq ($(MACHINE), baiz)
CC=pgcc
LAPACK=-llapack -lblas -pgf77libs    
LAPACKL=/opt/pgi/linux86-64/10.8/lib/
FFTW3= -lfftw3
FFTW3L= -L/usr/lib   
FFTW3I=    
CC_FLAGS= -O3 -lm
endif

### Compiler settings medusa
ifeq ($(MACHINE), medusa)
CC=/usr/local/intel/cc/9.0/bin/icc
LAPACK=-lmkl_lapack -lmkl_ia32 -lguide -lpthread
LAPACKL=-L/usr/local/intel/mkl/8.0.1/lib/32
FFTW3=-lfftw3
FFTW3L=-L/home/lacour/fftw3/lib
FFTW3I=-I/home/lacour/fftw3/include
CC_FLAGS= -O3 -lm
endif

### Compiler settings new medusa
ifeq ($(MACHINE), newmedusa)
CC=icc
LAPACK=-lmkl_intel_lp64 -lmkl_lapack -lmkl_core -lguide -lmkl_intel_thread -lpthread 
LAPACKL=-L/usr/local/intel/mkl10.2/lib/em64t
FFTW3=-lfftw3f
FFTW3L=-L/home/lacour/fftw312/lib
FFTW3I=-I/home/lacour/fftw312/include
CC_FLAGS= -O3 -lm
endif

### Compiler Settings Hydron
ifeq ($(MACHINE), hydron)
CC = pgcc
LAPACK_DIR=/usr/local/pgi/linux86-64/7.0-7/lib/
CC_FLAGS = -O3 -L/usr/local/lib -lfftw3 -lm
LAPACK = -llapack -lblas -pgf77libs
LAPACKL = -L$(LAPACK_DIR) 
endif

### Compiler Settings HPCIBM
ifeq ($(MACHINE), hpcibm)
CC = /opt/pgi/linux86/6.0/bin/pgcc
CPPFLAGS= -B -I/opt/fftw-3.1.1/include -L/opt/fftw-3.1.1/lib
CC_FLAGS = -B -O3 -L/opt/fftw-3.1.1/lib -lfftw3 -lm
LAPACK = -llapack -lblas -lpgftnrtl
endif

### Compiler Settings Millipede
ifeq ($(MACHINE), millipede)
CC = icc
CPPFLAGS= -B -I/cm/shared/apps/fftw/intel/64/3.2.2/include -L/cm/shared/apps/fftw/intel/64/3.2.2/lib
CC_FLAGS = -B -O3 -L/cm/shared/apps/fftw/intel/64/3.2.2/lib -lfftw3 -lm
LAPACK = -llapack -lblas -lifcore
LAPACKL = -L$(LAPACK_DIR) -L$(BLASDIR)
endif

### Compiler Settings Bohr
ifeq ($(MACHINE), bohr)
CC = icc
FFTW3=-lfftw3
FFTW3L=-L/software/apps/fftw/intel/3.3.3/lib -L/software/apps/intel/compiler/2013/13.0.1/compiler/lib/intel64 -L/software/apps/intel/compiler/2013/13.0.1/mkl/lib/intel64
FFTW3I=-I/software/apps/fftw/intel/3.3.3/include -openmp
CC_FLAGS = -B -O3 -lm -openmp
LAPACK = -mkl -lmkl_intel_lp64 -lmkl_core -liomp5 -lmkl_intel_thread -lpthread
#LAPACKL = -L/usr/lib64/atlas
#LAPACKL = -L$(LAPACK_DIR) -L$(BLASDIR)
endif


CPPFLAGS= $(FFTW3I)
LIBS=$(LAPACKL) $(LAPACK) $(FFTW3L) $(FFTW3)
#CC_FLAGS=-B -O3 -lm $(LAPACKL) $(LAPACK) $(FFTW3L) $(FFTW3)

all: echo NISE3 translate 2DFFT

echo:
	@echo Compiling on the machine $(MACHINE) with $(CC)

NISE3: NISE3.o NISE3subs.o
	$(CC) NISE3.o NISE3subs.o $(LIBS) -o NISE3 $(CC_FLAGS)

NISE3.o: NISE3.c NISE3.h typesA.h lapack.h

NISE3subs.o: NISE3subs.c NISE3subs.h typesA.h lapack.h

translate.o: translate.c translate.h

translate: translate.o nonsubs.o
	$(CC) translate.o nonsubs.o $(LIBS) -o translate $(CC_FLAGS)


nonsubs.o: nonsubs.c nonsubs.h types.h lapack.h

2DFFT: 2DFFT.o nrutil.o
	$(CC) nrutil.o 2DFFT.o $(LIBS) -o 2DFFT $(CC_FLAGS) 

2DFFT.o: 2DFFT.c

nrutil.o: nrutil.h nrutil.c

clean:
	\rm *.o
