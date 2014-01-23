SOURCES    = acor.cpp arctan.cpp covariance.cpp data.cpp ensemblemean.cpp exception.cpp exoplanet_hyperpara.cpp exoplanet_init.cpp exoplanetjd.cpp factorial.cpp int2str.cpp linearalgebra.cpp keplers_eqn.cpp mean.cpp model.cpp numzero.cpp orbcrocor.cpp randn.cpp record.cpp rng.cpp rosenbrock.cpp sampling.cpp samplingmc.cpp var.cpp

OBJECTS    = $(SOURCES:.cpp=.o) main.o

EXECUTABLE = MCMC

CC         = g++

TXT        = 

ALL        = $(SOURCES:.cpp=.h) $(SOURCES) Makefile main.cpp  $(TXT)

#CFLAGS = -g #for debugging
CFLAGS = -O3 #for running

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(CFLAGS) -o $(EXECUTABLE) $(OBJECTS)  -llapack -lblas

%.o: %.cpp
	$(CC) $(CFLAGS) -c $<

clean:
	rm $(EXECUTABLE) $(OBJECTS)
	
tarball: $(ALL)
	tar -cvf BIC.tar $(ALL)
