FC = pgf90

FLGS = -O2 -pc 64 -Mdalign -Mextend -Mnoframe -byteswapio -Bstatic_pgi

TARGT = GC_EXT_prep

OBJCTS = read_chem.o \
	main.o

$(TARGT) : $(OBJCTS)
	$(FC) $(FLGS) -o $(TARGT) $(OBJCTS)

.F.o :
	$(FC) -c -o $@ $(FLGS) $<

clean :
	rm *.o

