F2PY = f2py

SRC = $(wildcard *.f90)
# OBJ = $(SRC:.f90=)
OBJ = $(SRC:.f90=.so)

all: $(OBJ)
%.so: %.f90
	$(F2PY) -c  -m $(subst .so,,$@) $<
clean:
	rm $(OBJ)
