all: RE-inverse7 RE-cube7 RE-aes

RE-inverse7: ReductionExhaustive.cpp
	g++ -std=c++2a -DSBOX_INVERSE7 -O3 ReductionExhaustive.cpp -fopenmp -o RE-inverse7

RE-cube7: ReductionExhaustive.cpp
	g++ -std=c++2a -DSBOX_CUBE7 -O3 ReductionExhaustive.cpp -fopenmp -o RE-cube7

RE-aes: ReductionExhaustive.cpp
	g++ -std=c++2a -DSBOX_AES -O3 ReductionExhaustive.cpp -fopenmp -o RE-aes

RE-inverse8: ReductionExhaustive.cpp
	g++ -std=c++2a -DSBOX_INVERSE8 -O3 ReductionExhaustive.cpp -fopenmp -o RE-inverse8

RE-cube8: ReductionExhaustive.cpp
	g++ -std=c++2a -DSBOX_CUBE8 -O3 ReductionExhaustive.cpp -fopenmp -o RE-cube8

RE-apn8c1: ReductionExhaustive.cpp
	g++ -std=c++2a -DSBOX_APN8_C1 -O3 ReductionExhaustive.cpp -fopenmp -o RE-apn8c1

RE-five8: ReductionExhaustive.cpp
	g++ -std=c++2a -DSBOX_FIVE8 -O3 ReductionExhaustive.cpp -fopenmp -o RE-five8

RE-five9: ReductionExhaustive.cpp
	g++ -std=c++2a -DSBOX_FIVE9 -O3 ReductionExhaustive.cpp -fopenmp -o RE-five9

RE-random8: ReductionExhaustive.cpp
	g++ -std=c++2a -DSBOX_RANDOM8 -O3 ReductionExhaustive.cpp -fopenmp -o RE-random8

RE-custom: ReductionExhaustive.cpp
	g++ -std=c++2a -DSBOX_CUSTOM -O3 ReductionExhaustive.cpp -fopenmp -o RE-custom

RE-custom-large: ReductionExhaustive.cpp
	g++ -std=c++2a -DSBOX_CUSTOM -DSBOX_CUSTOM_LARGE -O3 ReductionExhaustive.cpp -fopenmp -o RE-custom-large