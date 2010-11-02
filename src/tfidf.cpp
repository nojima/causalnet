#include <cstdio>
#include <cmath>
#include <map>
#include <vector>
#include <algorithm>
using namespace std;

struct TfIdf {
  typedef map<int, int> Map;
  typedef Map::iterator Iter;

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
  Iter it = tf[docId].lower_bound(wordId);
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

  int nnz = 0;
  for (int docId = 0; docId < docCount; ++docId) {
    for (Iter it = tf[docId].begin(); it != tf[docId].end(); ) {
      if (minDf <= df[it->first] && df[it->first] <= maxDf) {
        ++nnz;
        ++it;
      } else {
        tf[docId].erase(it++);
      }
    }
  }

  fprintf(file, "%d\t%d\t%d\n\n", docCount, wordCount, nnz);

  for (int docId = 0; docId < docCount; ++docId) {
    double sum = 0.0;
    for (Iter it = tf[docId].begin(); it != tf[docId].end(); ++it) {
      double v = it->second * log((double)docCount / df[it->first]);
      sum += v * v;
    }
    double inv = 1.0 / sqrt(sum);
    for (Iter it = tf[docId].begin(); it != tf[docId].end(); ++it) {
      fprintf(file, "%.8f\n", inv * it->second * log((double)docCount / df[it->first]));
    }
  }
  fprintf(file, "\n");

  int* ptr = new int[docCount+1];
  ptr[0] = 0;
  for (int docId = 0; docId < docCount; ++docId) {
    int count = 0;
    for (Iter it = tf[docId].begin(); it != tf[docId].end(); ++it) {
      fprintf(file, "%d\n", it->first);
      ++count;
    }
    ptr[docId+1] = ptr[docId] + count;
  }
  fprintf(file, "\n");

  for (int docId = 0; docId < docCount + 1; ++docId) {
    fprintf(file, "%d\n", ptr[docId]);
  }
  fprintf(file, "\n");

  delete ptr;

  return 0;
}

int main()
{
  TfIdf tfidf;
  tfidf.load(stdin);
  tfidf.save(stdout, 5, tfidf.docCount * 0.1);
}
