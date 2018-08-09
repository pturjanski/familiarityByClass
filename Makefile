# Makefile for familiarityByClass

# options
CC		= g++
XXFLAGS		= -O3 -fomit-frame-pointer 
CPPFLAGS	= -Wall -DNDEBUG 
LDFLAGS		= 
LDLIBS		= -lboost_system -lboost_filesystem

# targets
.PHONY: all
all: familiarityByClass
familiarityByClass: sais.o Node.o Trie.o FastaElement.o Tools.o mr.o smrUnn.o familiarityByClass.o

distclean: clean
clean:
	$(RM) familiarityByClass familiartyChart.o mr.o smrUnn.o FastaElement.o Tools.o Trie.o Node.o sais.o 

# dependencies
sais.o Node.o Trie.o FastaElement.o Tools.o mr.o smrUnn.o familiarityByClass.o: sais.h Node.h Trie.h FastaElement.h Tools.h mr.h smrUnn.h Makefile
