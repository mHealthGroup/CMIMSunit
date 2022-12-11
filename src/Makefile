HEADERS = mims_aggregate.h stats/approx.h mims_combine_axes.h mims_extrapolate.h mims_filter.h \
	mims_helper.h stats/pgamma.h signal/butter.h signal/arma.h signal/bilinear.h \
	signal/poly.h signal/sftrans.h signal/filter.h stats/spline.h \
	stats/smspline.h trapz/trapz.h mims_unit.h

OBJECTS = mims_aggregate.o stats/approx.o mims_combine_axes.o mims_extrapolate.o mims_filter.o \
	mims_helper.o stats/pgamma.o signal/butter.o signal/arma.o signal/bilinear.o \
	signal/poly.o signal/sftrans.o signal/filter.o stats/spline.o \
	stats/smspline.o trapz/trapz.o mims_unit.o

default: mims_unit

%.o: %.c $(HEADERS)
	gcc -c -g $< -o $@

mims_unit: $(OBJECTS)
	gcc $(OBJECTS) -o $@

clean:
	-rm -f $(OBJECTS)
	-rm -f mims_unit
	-make