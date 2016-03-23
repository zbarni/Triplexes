CC=g++
CFLAGS=-c -Wall -g -O0

SOURCES=myers-main.cpp myers/edlib.cpp
LDFLAGS=
INCLUDES=-I/home/barni/code/rwth/bioinformatics/git/triplexator/ 
#INCLUDES=-I/home/zbarni/code/triplexator-project/seqan/seqan-library-1.4.0/include/
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=m

all: $(SOURCES) $(EXECUTABLE) 

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDES) $< -o $@

clean:
	rm *.o 
