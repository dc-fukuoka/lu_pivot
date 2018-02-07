fc      = ifort
fcflags = -D_DEBUG -g -O3 -mavx -fopenmp
srcs    = lu_pivot.F90
objs    = $(srcs:%.F90=%.o)
bin     = a.out
vsl     = mkl_vsl.mod
libs    = -mkl

all: $(bin)

$(vsl):  $(MKLROOT)/include/mkl_vsl.f90
	$(fc) -c $(fcflags) $<

$(objs): $(srcs) $(vsl)
	$(fc) -c $(fcflags) $<

$(bin): $(objs)
	$(fc) $(fcflags) $^ -o $@ $(libs)

clean:
	rm -f $(bin) *.mod
