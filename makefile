CFLAGS = -O3 -ffast-math -IDIR=/usr/include/g++-3/
LDFLAGS = -lstdc++
vpath %.o objects
vpath %.cpp source
vpath %.h source
objdir = objects
targets = SAXS_calc cylgen

all : $(targets)


%.o : %.cpp
	gcc -c $(CFLAGS) $(CPPFLAGS) $< -o objects/$@

SAXS_calc : SAXS_calc.o ScatteringObject.o PeriodicCylinderScatteringObject.o \
	PeriodicCylinderCustomAxialDensityScatteringObject.o \
	SpecialFunctionGenerator.o
	gcc $(FLAGS) $(LDFLAGS) -o SAXS_calc $(objdir)/SAXS_calc.o \
	$(objdir)/ScatteringObject.o $(objdir)/PeriodicCylinderScatteringObject.o \
	$(objdir)/PeriodicCylinderCustomAxialDensityScatteringObject.o \
	$(objdir)/SpecialFunctionGenerator.o


cylgen : CylinderGenerator.o 
	gcc $(FLAGS) $(LDFLAGS) -o cylgen $(objdir)/CylinderGenerator.o 



.PHONY : clean
clean :
	rm objects/*.o


.PHONY : delete
delete : 
	rm $(targets)






