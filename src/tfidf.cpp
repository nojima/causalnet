#include <cstdio>
#include <cmath>
#include <map>
#include <vector>
#include <algorithm>
using namespace std;

struct TfIdf {
  typedef map<int, int> Map;

  TfIdf(): docCount(0), wordCount(0), tf(NULL), df(NULL) {}
  ~TfIdf();
  void load(FILE* file);
  int save(FILE* file, int minDf, int maxDf) const;
  void addWord(int docId, int wordId);

  int docCount;       // the number of documents
  int wordCount;      // the number of words
  Map* tf;            // term frequency
  int* df;            // document frequency
};

TfIdf::~TfIdf()
{
  delete [] tf;
  delete [] df;
}

void TfIdf::load(FILE* file)
{
  fscanf(file, "%d%d", &docCount, &wordCount);
  tf = new Map[docCount];
  df = new int[wordCount];
  int docId, wordId;
  while (fscanf(file, "%d%d", &docId, &wordId) != EOF) {
    addWord(docId - 1, wordId - 1);
  }
}

void TfIdf::addWord(int docId, int wordId)
{
  Map::iterator it = tf[docId].lower_bound(wordId);
  if (it != tf[docId].end() && it->first == wordId) {
    it->second += 1;
  } else {
    tf[docId].insert(it, make_pair(wordId, 1));
    df[wordId] += 1;
  }
}

int TfIdf::save(FILE* file, int minDf, int maxDf) const
{
  if (file == NULL) return -1;
  fprintf(file, "%d\t%d\n", docCount, wordCount);
  vector<double> tfidf;
  vector<int> index;
  for (int docId = 0; docId < docCount; ++docId) {
    double sum = 0.0;
    index.clear();
    tfidf.clear();
    for (Map::iterator it = tf[docId].begin(); it != tf[docId].end(); ++it) {
      if (minDf <= df[it->first] && df[it->first] <= maxDf) {
        index.push_back(it->first);
        tfidf.push_back(it->second * log((double)docCount / df[it->first]));
        sum += tfidf.back();
      }
    }
    double inv = 1.0 / sum;
    for (size_t k = 0; k < tfidf.size(); ++k) {
      fprintf(file, "%d\t%d\t%f\n", docId, index[k], inv * tfidf[k]);
    }
  }
  return 0;
}

int main()
{
  TfIdf tfidf;
  tfidf.load(stdin);
  tfidf.save(stdout, 5, tfidf.docCount * 0.1);
}
