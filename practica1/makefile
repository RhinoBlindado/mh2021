INC = -I ./libs/

COPKM_exe : COPKM.o random.o
	g++ -std=c++17 $(INC) -O2 -o $@ $^
	
COPKM.o : COPKM.cpp
	g++ -std=c++17 $(INC) -O2 -c $< 

#BL_exe : BL.cpp random.o
#	g++ -std=c++17 -O2 -o $@ $<

random.o : libs/random.cpp
	g++ -std=c++17 $(INC) -c $<
