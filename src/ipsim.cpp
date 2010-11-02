#include <cstdio>
#include <map>
using namespace std;

class InnerProductSimilarity {
  public:
    InnerProductSimilarity(): n(0), m(0), a(NULL), sim(NULL) {}
    void load(FILE* file);
    void calc();
    int save(FILE* file) const;

  private:
    typedef map<int, double> Map;
    int n, m;
    Map* a;
    Map* sim;
};

void InnerProductSimilarity::load(FILE* file)
{
  fscanf(file, "%d%d", &n, &m);
  a = new Map[m];
  sim = new Map[n];
  int i, j;
  double aij;
  while (fscanf(file, "%d%d%lf", &i, &j, &aij) != EOF) {
    a[j][i] = aij;
  }
}

template<typename Iter>
inline Iter next(Iter it) { return ++it; }

void InnerProductSimilarity::calc()
{
  int nnz = 0;
  for (int j = 0; j < m; ++j) {
    for (Map::iterator it1 = a[j].begin(); it1 != a[j].end(); ++it1) {
      int i = it1->first;
      double aij = it1->second;
      for (Map::iterator it2 = next(it1); it2 != a[j].end(); ++it2) {
        int k = it2->first;
        double akj = it2->second;
        Map::iterator it3 = sim[i].lower_bound(k);
        if (it3 != sim[i].end() && it3->first == k) {
          it3->second += aij * akj;
        } else {
          sim[i].insert(it3, make_pair(k, aij * akj));
          ++nnz;
        }
      }
    }
    printf("nnz = %d\n", nnz);
  }
}

int InnerProductSimilarity::save(FILE* file) const
{
  if (file == NULL) { return -1; }
  fprintf(file, "%d\n", n);
  for (int i = 0; i < n; ++i) {
    for (Map::iterator it = sim[i].begin(); it != sim[i].end(); ++it) {
      fprintf(file, "%d\t%d\t%f\n", i, it->first, it->second);
      fprintf(file, "%d\t%d\t%f\n", it->first, i, it->second);
    }
  }
  return 0;
}

int main()
{
  InnerProductSimilarity ipsim;
  ipsim.load(stdin);
  ipsim.calc();
  ipsim.save(stdout);
}
