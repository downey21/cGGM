
# -*- coding: utf-8 -*-

rm(list = ls())

# install.packages("MRCE")

library(MRCE)

set.seed(48105)
n <- 50
p <- 10
q <- 5

Omega.inv <- diag(q)
for (i in 1:q) {
    for (j in 1:q) {
        Omega.inv[i,j] <- 0.7^abs(i-j)
    }
}
out <- eigen(Omega.inv, symmetric = TRUE)
Omega.inv.sqrt <- tcrossprod(out$vec*rep(out$val^(0.5), each = q), out$vec)
Omega <- tcrossprod(out$vec*rep(out$val^(-1), each = q), out$vec)

X <- matrix(rnorm(n*p), nrow = n, ncol = p)
E <- matrix(rnorm(n*q), nrow = n, ncol = q) %*% Omega.inv.sqrt
Beta <- matrix(stats::rbinom(p*q, size=1, prob=0.1)*runif(p*q, min=1, max=2), nrow = p, ncol = q)
mu <- 1:q

Y <- rep(1,n) %*% t(mu) + X %*% Beta + E

lam1.vec <- rev(10^seq(from = -2, to = 0, by = 0.5))
lam2.vec <- rev(10^seq(from = -2, to = 0, by = 0.5))
cvfit <- MRCE::mrce(Y = Y, X = X, lam1.vec = lam1.vec, lam2.vec = lam2.vec, method = "cv")
cvfit

fit <- MRCE::mrce(Y = Y, X = X, lam1 = 10^(-1.5), lam2 = 10^(-0.5), method = "single")
fit

lam2.mat <- 1000*(fit$Bhat==0)
refit <- MRCE::mrce(Y = Y, X = X, lam2 = lam2.mat, method = "fixed.omega", omega = fit$omega, tol.in = 1e-12) 
refit

fit$Bhat
refit$Bhat
