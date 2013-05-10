F2PY = f2py
FF = gfortran
FFLAGS = -shared -fPIC -O3

SRC = weno.f90 rp_sw1d_roe.f90
OBJ = $(SRC:.f90=.so)
FOBJ = $(SRC:.f90=_f2py.so)

SOLVERS = hr_solver1d.so adv_weno.so weno_solver1d.so hr_solver2d.so

all: $(SOLVERS) $(FOBJ) run

$(SOLVERS) : %.so : $(OBJ) %.f90 
	$(F2PY) -c -m $(subst .so,,$@) $^

$(FOBJ): %_f2py.so : %.f90
	$(F2PY) -c  -m $(subst .so,,$@) $<

$(OBJ) : %.so: %.f90
	$(FF) $(FFLAGS) -o $@ $<

run: hr_solver2d.so swe2d-lw.py
	python swe2d-lw.py
clean:
	rm $(OBJ) $(SRC:.f90=.mod) $(SOLVERS) $(FOBJ) 
