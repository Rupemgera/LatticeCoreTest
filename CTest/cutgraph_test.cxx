#include "LatticeCore/Core/parameterization/cut_graph.h"
#include "conf.h"
#include "gtest/gtest.h"
using namespace std;

TEST(cutgraph_test, test_spanning_tree)
{
  size_t level = 3;
  size_t m = 3;
  size_t n = (pow(m, level) - 1) / (m - 1);
  vector<LatticeCore::TreeEdge> edges;
  int u = 0;
  int k = 0;
  for (size_t i = 0; i * m < n; ++i) {
    for (size_t j = 1; j <= m && (i * m + j) < n; ++j) {
      edges.push_back({ i, i * m + j, k, k });
      k++;
    }
  }

  std::vector<std::pair<size_t, size_t>> query = {
    { 1, 6 }, { 4, 11 }, { 5, 2 }, { 2, 10 }, { 0, 12 }, { 0, 7 }, { 9, 7 }
  };

  LatticeCore::SpanningTree stree(n, edges);
  auto ans = stree.least_common_ancestor(query);
  for (int i = 0; i < ans.size(); ++i) {
    cout << query[i].first << "->" << query[i].second << " : ";
    cout << ans[i] << endl;
  }
  auto p = stree.path(4, 6, 1);
  for (auto& i : p) {
    cout << i << ' ';
  }
  cout << endl;
  auto order = stree.post_order();
  for (auto& i : order) {
    cout << i << ' ';
  }
  cout << endl;
}