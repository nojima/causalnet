CFLAGS = -Wall -Wno-unused-result -O2

all: tfidf ipsim
	
tfidf: src/tfidf.cpp
	g++ $(CFLAGS) src/tfidf.cpp -o tfidf

ipsim: src/ipsim.cpp
	g++ $(CFLAGS) src/ipsim.cpp -o ipsim

clean:
	rm tfidf ipsim
