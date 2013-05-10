F2PY = f2py
FF = gfortran
FFLAGS = -shared -fPIC 

SRC = weno.f90 rp_sw1d_roe.f90
# OBJ = $(SRC:.f90=)
OBJ = $(SRC:.f90=.so)
FOBJ = $(SRC:.f90=_f2py.so)
MODULE = hr_solver1d

all: ${MODULE} weno_solver1d.so adv_weno.so

adv_weno.so : $(OBJ) adv_weno.f90
	$(F2PY) -c -m adv_weno $^


weno_solver1d.so : $(OBJ) weno_solver1d.f90
	$(F2PY) -c -m weno_solver1d $^

${MODULE} : ${MODULE}.so

${MODULE}.so: $(OBJ) ${MODULE}.f90
	$(F2PY) -c -m ${MODULE}  $(OBJ) hr_solver1d.f90

%.so: %.f90
	$(F2PY) -c  -m $(subst .so,_f2py,$@) $<
	$(FF) $(FFLAGS) -o $@ $<

clean:
	rm $(OBJ) *.mod ${MODULE}.so $(FOBJ) weno_solver1d.so
