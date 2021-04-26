INC = -I ./libs/

all: GENE_exe MEME_exe

GENE_exe : GENE.o random.o
	g++ -std=c++17 $(INC) -O2 -o $@ $^
	
GENE.o : GENE.cpp
	g++ -std=c++17 $(INC) -O2 -c $< 

MEME_exe : MEME.o random.o
	g++ -std=c++17 $(INC) -O2 -o $@ $^
	
MEME.o : MEME.cpp
	g++ -std=c++17 $(INC) -O2 -c $< 

random.o : libs/random.cpp
	g++ -std=c++17 $(INC) -c $<

clean_obj:
	rm -f *.o

clean:
	rm -f *.o *_exe
