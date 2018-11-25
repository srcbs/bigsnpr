// [[Rcpp::depends(RcppParallel, BH, bigstatsr)]]
#include <bigstatsr/BMCodeAcc.h>
#include <RcppParallel.h>
using namespace RcppParallel;

struct Sum : public Worker {

  SubBMCode256Acc macc;
  double xySum;
  std::size_t j0, j;

  // constructors
  Sum(SubBMCode256Acc macc) :
    macc(macc),     xySum(0), j0(0), j(0) {}
  Sum(const Sum& sum, std::size_t j0, std::size_t j) :
    macc(sum.macc), xySum(0), j0(j0), j(j) {}
  Sum(const Sum& sum, Split) :
    macc(sum.macc), xySum(0), j0(sum.j0), j(sum.j) {}

  // accumulate just the element of the range I've been asked to
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; i++) {
      xySum += macc(i, j) * macc(i, j0);
    }
  }
  // join results
  void join(const Sum& rhs) {
    xySum += rhs.xySum;
  }
};

// [[Rcpp::export]]
NumericVector parallelVectorSum(Environment BM) {

  XPtr<FBM> xpBM = BM["address"];
  std::size_t n = xpBM->nrow();
  std::size_t m = xpBM->ncol();
  SubBMCode256Acc macc(xpBM, seq_len(n) - 1, seq_len(m) - 1, BM["code256"]);

  int grain = std::sqrt(n);

  Sum sum0(macc);
  NumericVector res(m);
  for (size_t j = 0; j < m; j++) {
    Sum sum(sum0, 0, j);
    parallelReduce(0, n, sum, grain);
    res[j] = sum.xySum;
  }

  return res;
}




/*** R
RcppParallel::setThreadOptions(2)
library(bigsnpr)
snp <- snp_attachExtdata()
G <- snp$genotypes
test0 <- parallelVectorSum(G)

G2 <- big_copy(G, ind.row = rep(rows_along(G), 500))
dim(G2)
RcppParallel::setThreadOptions(1)
system.time(test1 <- parallelVectorSum(G2))  # 100 / 3
testthat::expect_identical(test1, 500 * test0)
RcppParallel::setThreadOptions(2)
system.time(test2 <- parallelVectorSum(G2))  # 177 / 39
testthat::expect_identical(test2, 500 * test0)
# RcppParallel::setThreadOptions(4)
# system.time(parallelVectorSum(G2))  # 305 / 44
# RcppParallel::setThreadOptions(10)
# system.time(parallelVectorSum(G2))  # 689 / 53
*/
