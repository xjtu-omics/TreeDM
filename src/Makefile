
# TreeDM Makefile
#

CC		= g++
WARNINGS	= -Wall -ansi
DEBUG		= -ggdb

#OPTIMIZE       = -O3  -mcpu=x86-64 -mmmx -momit-leaf-frame-pointer -fomit-frame-pointer -ftracer -malign-double
#OPTIMIZE       = -O3 -march=pentium3 -msse -momit-leaf-frame-pointer -fomit-frame-pointer -ftracer -fprefetch-loop-arrays -malign-double
OPTIMIZE       = -O3 -march=x86-64 -msse -momit-leaf-frame-pointer -fomit-frame-pointer -fprefetch-loop-arrays -malign-double
# Athlon-XP
#OPTIMIZE        = -O3 -march=x86-64 -msse -mmmx -m3dnow -mfpmath=sse -momit-leaf-frame-pointer -fomit-frame-pointer -ftracer -fprefetch-loop-arrays -malign-double
#float optimized
#OPTIMIZE       = -O3 -march=x86-64  -momit-leaf-frame-pointer -fomit-frame-pointer -ftracer -fprefetch-loop-arrays -malign-double
#OPTIMIZE       = -m128bit-long-double -msse
PROFILE         = -pg
#CFLAGS         = $(WARNINGS)
CFLAGS          = $(WARNINGS) $(OPTIMIZE)
#CFLAGS         = $(WARNINGS) $(DEBUG)
#CFLAGS         = $(WARNINGS) $(PROFILE)
SHELL		= /bin/sh

IMPORTANTFILES	= Makefile PGOneI_PerProtein.cpp  PGOneI_PerProtein.h  PGOneI.cpp  PGOneI_PD.cpp  PGOneI_PD.h


all:	clean TreeDM

TreeDM: PGOneI_PerProtein.o PGOneI_PD.o PGOneI.o
	$(CC) $(CFLAGS) -o $@ $^



clean:
	rm -f *.o  TreeDM

backup:
	tar jvcf ./backup_TreeDM-`date +'%Y-%m-%d_%H.%M'`.tar.bz2 $(IMPORTANTFILES)




# Dependency's




PGOneI_PerProtein.o: PGOneI_PerProtein.cpp
	$(CC) $(CFLAGS) -c $^
	
PGOneI_PD.o: PGOneI_PD.cpp
	$(CC) $(CFLAGS) -c $^

PGOneI.o: PGOneI.cpp
	$(CC) $(CFLAGS) -c $^
	
