INSTALLDIR=/usr/local/bin

src = $(wildcard *.cpp)
hdr = $(wildcard *.h *.hpp lib/*.hpp)
libsrc = lib/delaunay.cpp lib/triangle.cpp

CC=g++
COPT=-g -std=c++11 -I. -I./lib
#COPT=-O3 --std=c++11

clipcli: $(src) $(hdr) $(libsrc)
	$(CC) $(COPT) -o $@ $^

.PHONY: clean
clean:
	rm -f *.o clipcli

.PHONY: install
install: clipcli
	mkdir -p $(INSTALLDIR)
	cp $< $(INSTALLDIR)

