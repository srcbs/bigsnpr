test <- snp_attachExtdata()
G <- test$genotypes
n <- nrow(G)
m <- ncol(G)

K <- 50
G2 <- big_copy(G, ind.col = rep(cols_along(G), K))

RcppParallel::setThreadOptions(1)
system.time(snp_clumping(G2, infos.chr = rep(1:K, each = m)))
system.time(snp_clumping(G2, infos.chr = rep(1:K, each = m), ncores = 2))
system.time(snp_clumping(G2, infos.chr = rep(1:K, each = m), ncores = 4))
system.time(snp_clumping(G2, infos.chr = rep(1:K, each = m), ncores = 10))

RcppParallel::setThreadOptions(2)
system.time(snp_clumping(G2, infos.chr = rep(1:K, each = m)))
RcppParallel::setThreadOptions(4)
system.time(snp_clumping(G2, infos.chr = rep(1:K, each = m)))
RcppParallel::setThreadOptions(10)
system.time(snp_clumping(G2, infos.chr = rep(1:K, each = m)))
# 65 // 36.5 // 21 // 14
# vs     35 // 20.5 // 12.5
