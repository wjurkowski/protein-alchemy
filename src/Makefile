#protein_alchemy 1.0 makefile szalu nima

F95 = gfortran-4.3 -ffree-form
OBJECTS=tablice.o main.o readparam.o readalign.o readseq.o readstr.o transform.o fikpsik.o winkiel.o dihsort.o ellipse.o structure.o eceppanal.o ramadist.o volume.o analyse_bs.o fitting.o fitting_capri.o grid.o sanal.o capri_rms.o vr.o crd2pdb.o blocker.o finanal.o functions.o

siatami: $(OBJECTS) 

	$(F95) -o pr_alchem $(OBJECTS) /usr/lib/liblapack.so.3gf
#	cp pr_alchem /usr/local/bin/

tablice.o: tablice.f95
	$(F95) -c tablice.f95
main.o: main.f95
	$(F95) -c main.f95
readparam.o: readparam.f95
	$(F95) -c readparam.f95
readalign.o: readalign.f95
	$(F95) -c readalign.f95
readseq.o: readseq.f95
	$(F95) -c readseq.f95	
readstr.o: readstr.f95
	$(F95) -c readstr.f95
transform.o: transform.f95
	$(F95) -c transform.f95
fikpsik.o: fikpsik.f95
	$(F95) -c fikpsik.f95 
winkiel.o: winkiel.f95
	$(F95) -c winkiel.f95
dihsort.o: dihsort.f95
	$(F95) -c dihsort.f95
ellipse.o: ellipse.f95
	$(F95) -c ellipse.f95
structure.o: structure.f95
	$(F95) -c structure.f95	
eceppanal.o: eceppanal.f95
	$(F95) -c eceppanal.f95
ramadist.o: ramadist.f95
	$(F95) -c ramadist.f95
volume.o: volume.f95
	$(F95) -c volume.f95
analyse_bs.o: analyse_bs.f95
	$(F95) -c analyse_bs.f95
fitting.o: fitting.f95
	$(F95) -c fitting.f95
fitting_capri.o: fitting_capri.f95
	$(F95) -c fitting_capri.f95
grid.o: grid.f95
	$(F95) -c grid.f95
sanal.o: sanal.f95
	$(F95) -c sanal.f95
capri_rms.o: capri_rms.f95
	$(F95) -c capri_rms.f95
vr.o: vr.f95
	$(F95) -c vr.f95	
crd2pdb.o: crd2pdb.f95
	$(F95) -c crd2pdb.f95
blocker.o: blocker.f95
	$(F95) -c blocker.f95
finanal.o: finanal.f95
	$(F95) -c finanal.f95
functions.o: functions.f95
	$(F95) -c functions.f95

clean: 	
	rm -f *.o *~ pr_alchem
