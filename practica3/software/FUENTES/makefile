INC = -I ./libs/

all: ES_exe BMB_exe ILS_exe

ES_exe : ES.o random.o 
	g++ -std=c++17 $(INC) -g -o $@ $^
	
ES.o : ES.cpp
	g++ -std=c++17 $(INC) -g -c $^ 

BMB_exe : BMB.o random.o
	g++ -std=c++17 $(INC) -O2 -o $@ $^
	
BMB.o : BMB.cpp
	g++ -std=c++17 $(INC) -O2 -c $< 

ILS_exe : ILS.o random.o
	g++ -std=c++17 $(INC) -O2 -o $@ $^
	
ILS.o : ILS.cpp
	g++ -std=c++17 $(INC) -O2 -c $< 	
	
random.o : libs/random.cpp
	g++ -std=c++17 $(INC) -c $<

clean_obj:
	rm -f *.o

clean:
	rm -f *.o *_exe
