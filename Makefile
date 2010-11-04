CFLAGS = -Wall -Wno-unused-result -O2

all: tfidf ipsim ap
	
tfidf: src/tfidf.cpp
	g++ $(CFLAGS) src/tfidf.cpp -o tfidf

ipsim: src/ipsim.cpp
	g++ $(CFLAGS) src/ipsim.cpp -o ipsim

ap: src/ap.cpp
	g++ $(CFLAGS) src/ap.cpp -o ap

clean:
	rm -f tfidf ipsim ap
