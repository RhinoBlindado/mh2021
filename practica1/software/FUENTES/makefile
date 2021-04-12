INC = -I ./libs/

all: COPKM_exe BL_exe

COPKM_exe : COPKM.o random.o
	g++ -std=c++17 $(INC) -O2 -o $@ $^
	
COPKM.o : COPKM.cpp
	g++ -std=c++17 $(INC) -O2 -c $< 

BL_exe : BL.o random.o
	g++ -std=c++17 $(INC) -O2 -o $@ $^
	
BL.o : BL.cpp
	g++ -std=c++17 $(INC) -O2 -c $< 

random.o : libs/random.cpp
	g++ -std=c++17 $(INC) -c $<

clean_obj:
	rm -f *.o

clean:
	rm -f *.o *_exe
