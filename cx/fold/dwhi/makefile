FCOMP = gfortran -c
FLINK = gfortran

OBJ = dwhi.o input.o pots.o elastic.o inelastic.o radint_4.o angmom.o coulomb.o subs.o


dwhi.o:		dwhi.for
	$(FCOMP) dwhi.for
input.o:	input.for
	$(FCOMP) input.for
pots.o:		pots.for
	$(FCOMP) pots.for
elastic.o:	elastic.for
	$(FCOMP) elastic.for
inelastic.o:	inelastic.for
	$(FCOMP) inelastic.for
radint_4.o:	radint_4.for
	$(FCOMP) radint_4.for
angmom.o:	angmom.for
	$(FCOMP) angmom.for
coulomb.o:	coulomb.for
	$(FCOMP) coulomb.for
subs.o:		subs.for
	$(FCOMP) subs.for
dwhi: dwhi.o input.o pots.o elastic.o inelastic.o radint_4.o angmom.o coulomb.o subs.o
	$(FLINK) $(OBJ) -o dwhi

