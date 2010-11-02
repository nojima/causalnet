#include <cstdio>
#include <map>
#include <vector>
using namespace std;

struct InnerProductSimilarity {
  void load(FILE* file);
  void calc(double minSim);
  int save(FILE* file) const;

  int docCount, wordCount, nnz;
  double* tfidf;
  int* wordInd;
  int* docPtr;

  vector<double> sim;
  vector<int> colInd;
  vector<int> rowPtr;
};

void InnerProductSimilarity::load(FILE* file)
{
  fscanf(file, "%d%d%d", &docCount, &wordCount, &nnz);
  tfidf = new double[nnz];
  wordInd = new int[nnz];
  docPtr = new int[docCount+1];
  for (int k = 0; k < nnz; ++k) { fscanf(file, "%lf", &tfidf[k]); }
  for (int k = 0; k < nnz; ++k) { fscanf(file, "%d", &wordInd[k]); }
  for (int k = 0; k <= docCount; ++k) { fscanf(file, "%d", &docPtr[k]); }
}

void InnerProductSimilarity::calc(double minSim)
{
  // transpose
  double* tfidf2 = new double[nnz];
  int* docInd = new int[nnz];
  int* wordPtr = new int[wordCount+1];
  int* freq = new int[wordCount];
  fill(freq, freq + wordCount, 0);
  for (int i = 0; i < nnz; ++i) { freq[wordInd[i]] += 1; }
  wordPtr[0] = 0;
  for (int wordId = 0; wordId < wordCount; ++wordId) {
    wordPtr[wordId+1] = wordPtr[wordId] + freq[wordId];
  }
  fill(freq, freq + wordCount, 0);
  for (int docId = 0; docId < docCount; ++docId) {
    for (int i = docPtr[docId]; i < docPtr[docId+1]; ++i) {
      int wordId = wordInd[i];
      int index = wordPtr[wordId] + freq[wordId];
      tfidf2[index] = tfidf[i];
      docInd[index] = docId;
      freq[wordId] += 1;
    }
  }
  delete [] freq;

  // inner product
  typedef map<int, double> Map;
  Map ip;
  rowPtr.push_back(0);
  for (int docId = 0; docId < docCount; ++docId) {
    if (docId > 0) fprintf(stderr, "\r");
    fprintf(stderr, "(%d%%) %d/%d, nnz = %zd", 100*(docId+1)/docCount, docId+1, docCount, sim.size());

    ip.clear();
    for (int i = docPtr[docId]; i < docPtr[docId+1]; ++i) {
      int wordId = wordInd[i];
      double v1 = tfidf[i];
      for (int j = wordPtr[wordId]; j < wordPtr[wordId+1]; ++j) {
        double v2 = tfidf2[j];
        ip[docInd[j]] += v1 * v2;
      }
    }
    int count = 0;
    for (Map::iterator it = ip.begin(); it != ip.end(); ++it) {
      if (it->second >= minSim) {
        colInd.push_back(it->first);
        sim.push_back(it->second);
        ++count;
      }
    }
    rowPtr.push_back(rowPtr.back() + count);
  }
  fprintf(stderr, "\n");

  delete [] tfidf2;
  delete [] docInd;
  delete [] wordPtr;
}

int InnerProductSimilarity::save(FILE* file) const
{
  if (file == NULL) { return -1; }
  fprintf(file, "%d %d\n\n", docCount, rowPtr[docCount]);
  for (int i = 0; i < rowPtr[docCount]; ++i) { fprintf(file, "%f\n", sim[i]); }
  fprintf(file, "\n");
  for (int i = 0; i < rowPtr[docCount]; ++i) { fprintf(file, "%d\n", colInd[i]); }
  fprintf(file, "\n");
  for (int i = 0; i <= docCount; ++i) { fprintf(file, "%d\n", rowPtr[i]); }
  fprintf(file, "\n");
  return 0;
}

int main()
{
  InnerProductSimilarity ipsim;
  ipsim.load(stdin);
  ipsim.calc(0.1);
  ipsim.save(stdout);
}
