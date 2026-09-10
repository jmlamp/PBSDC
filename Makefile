# Makefile for the proximal bundle method PBSDC

FF = gfortran 
#FF = gfortran -O2
FFLAGS = -fbounds-check -Wall
#OPEN = -fopenmp 
OPEN =
# -pg, -fopenmp 

all: tpbsdc

tpbsdc: constants.o bundle1.o bundle2.o functions.o norm_min.o pbsdc.o dbdc.o pbdc.o tpbsdc.o fun.o plqdf1.o pvmm.o pqsubs.o mqsubs.o
	$(FF) -o tpbsdc $(FFLAGS) $(OPEN) constants.o bundle1.o bundle2.o functions.o norm_min.o pbsdc.o dbdc.o pbdc.o tpbsdc.o fun.o plqdf1.o pvmm.o pqsubs.o mqsubs.o 

constants.mod: constants.o constants.f95
	$(FF) -c $(FFLAGS) $(OPEN) constants.f95
	
constants.o: constants.f95
	$(FF) -c $(FFLAGS) $(OPEN) constants.f95

bundle1.mod: constants.mod bundle1.o bundle1.f95 
	$(FF) -c $(FFLAGS) $(OPEN) bundle1.f95
	
bundle1.o: constants.mod bundle1.f95
	$(FF) -c $(FFLAGS) $(OPEN) bundle1.f95 

bundle2.mod: constants.mod bundle2.o bundle2.f95
	$(FF) -c $(FFLAGS) $(OPEN) bundle2.f95 
	
bundle2.o: constants.mod bundle2.f95
	$(FF) -c $(FFLAGS) $(OPEN) bundle2.f95 

functions.mod: constants.mod functions.o functions.f95
	$(FF) -c $(FFLAGS) $(OPEN) functions.f95 
	
functions.o: constants.mod functions.f95
	$(FF) -c $(FFLAGS) $(OPEN) functions.f95 

norm_min.mod: constants.mod bundle1.mod bundle2.mod functions.mod norm_min.o norm_min.f95
	$(FF) -c $(FFLAGS) $(OPEN) norm_min.f95 
	
norm_min.o: constants.mod bundle1.mod bundle2.mod functions.mod norm_min.f95
	$(FF) -c $(FFLAGS) $(OPEN) norm_min.f95 

pbsdc.mod: constants.mod bundle1.mod bundle2.mod functions.mod norm_min.mod pbsdc.o pbsdc.f95 
	$(FF) -c $(FFLAGS) $(OPEN) pbsdc.f95	 
	
pbsdc.o: constants.mod bundle1.mod bundle2.mod functions.mod norm_min.mod pbsdc.f95
	$(FF) -c $(FFLAGS) $(OPEN) pbsdc.f95 
	
dbdc.mod: constants.mod bundle1.mod bundle2.mod functions.mod dbdc.o dbdc.f95 
	$(FF) -c $(FFLAGS) $(OPEN) dbdc.f95	 
	
dbdc.o: constants.mod bundle1.mod bundle2.mod functions.mod dbdc.f95
	$(FF) -c $(FFLAGS) $(OPEN) dbdc.f95 
	
pbdc.mod: constants.mod bundle1.mod bundle2.mod functions.mod norm_min.mod pbdc.o pbdc.f95 
	$(FF) -c $(FFLAGS) $(OPEN) pbdc.f95 

pbdc.o: constants.mod bundle1.mod bundle2.mod functions.mod norm_min.mod pbdc.f95
	$(FF) -c $(FFLAGS) $(OPEN) pbdc.f95
	
tpbsdc.o: constants.mod bundle1.mod bundle2.mod functions.mod norm_min.mod pbsdc.mod tpbsdc.f95
	$(FF) -c $(FFLAGS) $(OPEN) tpbsdc.f95 

fun.o: constants.mod norm_min.mod fun.f95
	$(FF) -c $(FFLAGS) $(OPEN) fun.f95 

plqdf1.o: plqdf1.f
	$(FF) -c $(FFLAGS) $(OPEN) plqdf1.f 

pqsubs.o: pqsubs.f
	$(FF) -c $(FFLAGS) $(OPEN) pqsubs.f 

mqsubs.o: mqsubs.f
	$(FF) -c $(FFLAGS) $(OPEN) mqsubs.f 

pvmm.o: pvmm.f
	$(FF) -c $(FFLAGS) $(OPEN) pvmm.f 	

clean:	
	rm tpbsdc constants.mod constants.o bundle1.mod bundle1.o bundle2.mod bundle2.o functions.mod functions.o dbdc.mod dbdc.o pbdc.mod pbdc.o norm_min.mod norm_min.o pbsdc.mod pbsdc.o tpbsdc.o fun.o plqdf1.o pvmm.o pqsubs.o mqsubs.o 
	echo Clean done