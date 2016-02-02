CC=g++
CFLAGS=-c -Wall

SOURCES=main.cpp
LDFLAGS=
INCLUDES=-I/home/zbarni/code/triplexator-project/seqan/seqan-library-1.4.0/include/
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=qgram

all: $(SOURCES) $(EXECUTABLE) 

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDES) $< -o $@

clean:
	rm *.o 
