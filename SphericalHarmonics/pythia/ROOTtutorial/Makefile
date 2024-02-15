datfiles = $(filter-out line.dat,$(wildcard *.dat))

all: make_2D.so

make_2D.so:
	root -l -b -e '.L ./make_2D.C+' -q

clean:
	rm -f *_C.{d,so,dylib} *_C_ACLiC*.pcm

distclean: clean
ifneq ($(strip $(datfiles)),)
	rm -f $(datfiles)
endif
