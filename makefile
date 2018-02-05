fc      = ifort
fcflags = -D_DEBUG -g -O3 -mavx -fopenmp
src     = lu_pivot.F90
bin     = a.out

$(bin): $(src)
	$(fc) $(fcflags) $^ -o $@

all: $(bin)

clean:
	rm -f $(bin) *.mod
