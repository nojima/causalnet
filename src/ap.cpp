// An Implementation of Affinity Propergation
// See: Clustering by Passing Messages Between Data Points

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <algorithm>
#include <cassert>
using namespace std;

namespace {
  struct Edge {
    int src;      // index of source
    int dst;      // index of destination
    double s;     // similarity s(src, dst)
    double r;     // responsibility r(src, dst)
    double a;     // availability a(src, dst)

    bool operator<(const Edge& rhs) const { return s < rhs.s; }
  };

  typedef vector<Edge*> Edges;

  struct Graph {
    int n;                // the number of vertices
    int nnz;              // the number of edges
    Edges* outEdges;      // array of out edges of corresponding vertices
    Edges* inEdges;       // array of in edges of corresponding vertices
    Edge* edges;          // all edges
    Edge* preferences;    // all preferences
  };

  // Build graph from sparse similarity matrix stored in CSR format.
  // Note that this function does not check any errors in the given input.
  // Parameter:
  //   input: Input file handle.
  //   prefType:
  //     1: use median of similarities as preference
  //     2: use minimum of similarities as preference
  //     3: use min - (max - min) of similarities as preference
  Graph* buildGraph(FILE* input, int prefType)
  {
    int n, nnz;
    fscanf(input, "%d%d", &n, &nnz);
    Graph* graph = new Graph;
    graph->n = n;
    graph->nnz = nnz;
    graph->outEdges = new Edges[n];
    graph->inEdges = new Edges[n];
    graph->edges = new Edge[nnz];
    Edge* edges = graph->edges;
    memset(edges, 0, sizeof(Edge) * nnz);

    // read similarity matrix
    for (int r = 0; r < nnz; ++r) {
      fscanf(input, "%lf", &edges[r].s);
    }
    for (int r = 0; r < nnz; ++r) {
      fscanf(input, "%d", &edges[r].src);
    }
    int k = 0;
    for (int r = 0; r <= n; ++r) {
      int next;
      fscanf(input, "%d", &next);
      for (; k < next; ++k) {
        edges[k].dst = r - 1;
      }
    }

    // calculate preferences
    double pref;
    if (prefType == 1) {
      sort(edges, edges + nnz);
      pref = (nnz % 2) ? edges[nnz/2].s : (edges[nnz/2 - 1].s + edges[nnz/2].s) / 2.0;
    } else if (prefType == 2) {
      pref = min_element(edges, edges + nnz)->s;
    } else if (prefType == 3) {
      double minValue = min_element(edges, edges + nnz)->s;
      double maxValue = max_element(edges, edges + nnz)->s;
      pref = 2*minValue - maxValue;
    } else {
      assert(false);      // invalid prefType
    }
    graph->preferences = new Edge[n];
    Edge* preferences = graph->preferences;
    memset(preferences, 0, sizeof(Edge) * n);
    for (int i = 0; i < n; ++i) {
      preferences[i].s = pref;
      preferences[i].src = i;
      preferences[i].dst = i;
    }

    for (int i = 0; i < nnz + n; ++i) {
      Edge* p;
      if (i < nnz) {
        // similarity
        p = &edges[i];
        if (p->src == p->dst) { continue; }
      } else {
        // preference
        p = &preferences[i - nnz];
      }
      // add small noise to avoid degeneracies
      p->s += (1e-16 * p->s + 1e-300) * (rand() / (RAND_MAX + 1.0));
      // add out/in edges to vertices
      graph->outEdges[p->src].push_back(p);
      graph->inEdges[p->dst].push_back(p);
    }

    return graph;
  }

  void destroyGraph(Graph* graph)
  {
    delete [] graph->outEdges;
    delete [] graph->inEdges;
    delete [] graph->edges;
    delete [] graph->preferences;
    delete graph;
  }

  int ccc;
  double totalMessage;

  inline void update(double& variable, double newValue, double damping)
  {
    if (fabs(newValue - variable) > 0.1) ++ccc;
    totalMessage += fabs(newValue - variable);
    variable = damping * variable + (1.0 - damping) * newValue;
  }

  void updateResponsibilities(Graph* graph, double damping)
  {
    for (int i = 0; i < graph->n; ++i) {
      Edges& edges = graph->outEdges[i];
      int m = edges.size();
      double max1 = -HUGE_VAL, max2 = -HUGE_VAL;
      double argmax1 = -1;
      for (int k = 0; k < m; ++k) {
        double value = edges[k]->s + edges[k]->a;
        if (value > max1) { swap(max1, value); argmax1 = k; }
        if (value > max2) { max2 = value; }
      }
      // update responsibilities
      for (int k = 0; k < m; ++k) {
        if (k != argmax1) {
          update(edges[k]->r, edges[k]->s - max1, damping);
        } else {
          update(edges[k]->r, edges[k]->s - max2, damping);
        }
      }
    }
  }

  void updateAvailabilities(Graph* graph, double damping)
  {
    for (int k = 0; k < graph->n; ++k) {
      Edges& edges = graph->inEdges[k];
      int m = edges.size();
      // calculate sum of positive responsibilities
      double sum = 0.0;
      for (int i = 0; i < m-1; ++i) {
        sum += max(0.0, edges[i]->r);
      }
      // calculate availabilities
      double rkk = edges[m-1]->r;
      for (int i = 0; i < m-1; ++i) {
        update(edges[i]->a, min(0.0, rkk + sum - max(0.0, edges[i]->r)), damping);
      }
      // calculate self-availability
      update(edges[m-1]->a, sum, damping);
    }
  }

  int updateExamplars(Graph* graph, vector<int>& examplar)
  {
    int changed = 0;
    for (int i = 0; i < graph->n; ++i) {
      Edges& edges = graph->outEdges[i];
      int m = edges.size();
      double maxValue = -HUGE_VAL;
      int argmax = i;
      for (int k = 0; k < m; ++k) {
        double value = edges[k]->a + edges[k]->r;
        if (value > maxValue) {
          maxValue = value;
          argmax = edges[k]->dst;
        }
      }
      if (examplar[i] != argmax) {
        examplar[i] = argmax;
        ++changed;
      }
    }
    return changed;
  }
}

// Cluster data points with Affinity Propagation.
// Parameters:
//   input: Input file which contains sparse similarity matrix. see buildGraph().
//   prefType: Specify what kind of preference we use. see buildGraph().
//   damping: The damping factor. (0.5 <= damping < 1.0)
//   maxit: The maximum number of iterations.
//   convit: Specify how many iterations this algorithm stops when examplars
//           did not change for.
// Returns:
//   Array of examplars of corresponding data points.
vector<int> affinityPropagation(FILE* input, int prefType, double damping, int maxit, int convit)
{
  assert(0.499 < damping && damping < 1.0);

  fprintf(stderr, "loading...\n");
  Graph* graph = buildGraph(input, prefType);
  vector<int> examplar(graph->n, -1);

  for (int i = 0, nochange = 0; i < maxit && nochange < convit; ++i, ++nochange) {
    ccc = 0;
    totalMessage = 0.0;
    
    updateResponsibilities(graph, damping);
    updateAvailabilities(graph, damping);
    int changed = updateExamplars(graph, examplar);
    if (changed) { nochange = 0; }

    fprintf(stderr, "iteration #%d (changed = %d, total message = %f, ccc = %d)\n", i, changed, totalMessage, ccc);
  }
  
  destroyGraph(graph);
  return examplar;
}

int main()
{
  vector<int> examplar = affinityPropagation(stdin, 1, 0.95, 3000, 20);
  printf("%zd\n\n", examplar.size());
  for (size_t i = 0; i < examplar.size(); ++i) {
    printf("%d\n", examplar[i]);
  }
}
