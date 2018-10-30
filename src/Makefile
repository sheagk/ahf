#******************************************************************#
#                                                                  #
# Makefile for               AHF                                   #
#                           =====                                  #
#                                                                  #
# version 1.0 (30/01/2005), Alexander Knebe                        #
#                                                                  #
# AHF (AMIGA-Halo-Finder) added by Stuart Gill (30/01/2005)        #
# 
# Makefile replaced by Benjamin Kimock (2018-05-21)
#                                                                  #
#******************************************************************#

.PHONY: clean veryclean

#*******************************************************************#
DEFINEFLAGS = ${MASTER_DEFINEFLAGS} -DAHF
export DEFINEFLAGS # probably not needed anymore

CFILES := $(wildcard *.c **/*.c)
OFILES := $(CFILES:%.c=%.o)

CFLAGS	+= $(OPTIMIZE) $(CCFLAGS) $(DEFINEFLAGS)

AHF: $(OFILES)
	$(CC) $(CFLAGS) $(LNFLAGS) -o ../bin/AHF-v1.0-094 $(OFILES)

# Replace the default rule for object files
%.o : %.c
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f $(OFILES)
	rm ../bin/AHF-v1.0-094

