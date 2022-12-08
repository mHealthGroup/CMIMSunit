HEADERS = aggregate.h approx/stats_approx.h combine_axes.h extrapolate.h filter.h helper.h pgamma/stats_pgamma.h \
	signal/butter/signal_butter.h signal/butter/arma.h signal/butter/bilinear.h signal/butter/poly.h \
	signal/butter/sftrans.h signal/filter/signal_filter.h spline/stats_spline.h \
	spline/stats_smspline.h trapz/trapz.h mims_unit.h

OBJECTS = aggregate.o approx/stats_approx.o combine_axes.o extrapolate.o filter.o helper.o pgamma/stats_pgamma.o \
	signal/butter/signal_butter.o signal/butter/arma.o signal/butter/bilinear.o signal/butter/poly.o \
	signal/butter/sftrans.o signal/filter/signal_filter.o spline/stats_spline.o \
	spline/stats_smspline.o trapz/trapz.o mims_unit.o

default: mims_unit

%.o: %.c $(HEADERS)
	gcc -c -g $< -o $@

mims_unit: $(OBJECTS)
	gcc $(OBJECTS) -o $@

clean:
	-rm -f $(OBJECTS)
	-rm -f mims_unit
	-make