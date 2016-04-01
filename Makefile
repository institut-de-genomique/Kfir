BOOST = /usr/include/boost/

CLASS_FILE=DnaDictionary.cpp ReadFile.cpp
OBJ = ${CLASS_FILE:.cpp=.o}
OBJ += dust.o
MAIN_FILE=kfir.cpp
BIN=kfir

COMPILER=g++
CPPFLAGS += -I$(BOOST) -I. -I$(SEQAN)
CPPFLAGS += -O3
CPPFLAGS += -pedantic -W -Wall -Wno-long-long -Wno-variadic-macros
CPPFLAGS += -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64
CPPFLAGS += -g
LDFLAGS = -L. -lgzstream -lz
AR = ar cr

COMPILER_BIS=gcc
CFLAGS += -I.
CFLAGS += -O3
CFLAGS += -c


kfir: libgzstream.a $(OBJ)
	${COMPILER} ${CPPFLAGS} -o $(BIN) $(MAIN_FILE) $(OBJ) ${LDFLAGS}

gzstream.o:
	${COMPILER} ${CPPFLAGS} -c -o gzstream.o gzstream.cpp

libgzstream.a: gzstream.o
	${AR} libgzstream.a gzstream.o

dust.o:
	${COMPILER_BIS} ${CFLAGS} -o dust.o dust.c

clean:
	rm -rf *.o a.out

