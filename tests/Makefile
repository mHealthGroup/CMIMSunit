HEADERS = ../src/mims_aggregate.h ../src/stats/approx.h ../src/mims_combine_axes.h ../src/mims_extrapolate.h ../src/mims_filter.h \
	../src/mims_helper.h ../src/stats/pgamma.h ../src/signal/filter.h ../src/stats/spline.h \
	../src/stats/smspline.h ../src/trapz/trapz.h ../src/mims_unit.h basic_tests.h

OBJECTS = ../src/mims_aggregate.o ../src/stats/approx.o ../src/mims_combine_axes.o ../src/mims_extrapolate.o ../src/mims_filter.o \
	../src/mims_helper.o ../src/stats/pgamma.o ../src/signal/filter.o ../src/stats/spline.o \
	../src/stats/smspline.o ../src/trapz/trapz.o ../src/mims_unit.o basic_tests.o demo.o

default: demo

%.o: %.c $(HEADERS)
	gcc -c -g $< -o $@

demo: $(OBJECTS)
	gcc $(OBJECTS) -o $@

clean:
	-rm -f $(OBJECTS)
	-rm -f demo
	-make