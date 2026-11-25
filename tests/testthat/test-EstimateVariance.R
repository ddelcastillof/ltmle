context("Test variance estimation")

#checks robust sd (tmle, iptw) is >= ic sd, and larger for tmle, under positivity violation
ExpectRobustStdDevUnderPositivity <- function(MakeY) {
  niter <- 40
  std.dev.ratio.tmle <- std.dev.ratio.iptw <- numeric(niter)
  for (i in 1:niter) {
    n <- 1000
    W <- rnorm(n)
    A <- rbinom(n, 1, plogis(5*W))
    Y <- MakeY(n)
    expect_warning({
      r.var.tmle <- ltmle(data.frame(W, A, Y), Anodes="A", Ynodes="Y", abar = 1, estimate.time=FALSE)
      r.var.ic <- ltmle(data.frame(W, A, Y), Anodes="A", Ynodes="Y", abar = 1, estimate.time=FALSE, variance.method="ic")
      r.var.iptw <- ltmle(data.frame(W, A, Y), Anodes="A", Ynodes="Y", abar = 1, estimate.time=FALSE, variance.method="iptw")}
      , "Variance estimate is based on")
    std.dev.ratio.tmle[i] <- summary(r.var.tmle)$treatment$std.dev / summary(r.var.ic)$treatment$std.dev
    std.dev.ratio.iptw[i] <- summary(r.var.iptw)$treatment$std.dev / summary(r.var.ic)$treatment$std.dev
  }
  expect_true(mean(std.dev.ratio.tmle) > 1.2)
  expect_true(mean(std.dev.ratio.tmle) > mean(std.dev.ratio.iptw))
  expect_true(min(std.dev.ratio.iptw) >= 1)
  expect_true(min(std.dev.ratio.tmle) >= 1)
}

test_that("TMLE based variance estimate > IC based variance estimate under positivity",{
  skip_on_cran()
  ExpectRobustStdDevUnderPositivity(function(n) rbinom(n, 1, 0.5))
})

test_that("All variance estimation methods run", {
  n <- 300
  W1 <- rnorm(n)
  W2 <- rnorm(n)
  A1 <- rexpit(W1 + W2)
  L1 <- rnorm(n) + A1
  C1 <- rexpit(W1 + W2)
  Y1 <- rexpit(W1 + W2 + A1)
  A2 <- rexpit(W1 + W2 + A1)
  L2 <- rnorm(n) + A1 + L1
  C2 <- rexpit(W1 + W2 + L1 + L2)
  Y2 <- rexpit(W1 + W2 - L1 - L2 + A1 - 2*A2)
  Y2[Y1 == 1] <- 1
  data <- data.frame(W1, W2, A1, L1, C1=BinaryToCensoring(is.uncensored = C1), Y1, A2, L2, C2=BinaryToCensoring(is.uncensored = C2), Y2)
  Anodes <- c("A1", "A2")
  Lnodes <- c("L1", "L2")
  Cnodes <- c("C1", "C2")
  Ynodes <- c("Y1", "Y2")
  r1 <- ltmle(data=data, Anodes=Anodes, Lnodes=Lnodes, Cnodes=Cnodes, Ynodes=Ynodes, abar=c(0,1), estimate.time = F, survivalOutcome = T)
  r2 <- ltmle(data=data, Anodes=Anodes, Lnodes=Lnodes, Cnodes=Cnodes, Ynodes=Ynodes, abar=c(0,1), estimate.time = F, survivalOutcome = T, variance.method = "iptw")
  expect_warning(r3 <- ltmle(data=data, Anodes=Anodes, Lnodes=Lnodes, Cnodes=Cnodes, Ynodes=Ynodes, abar=c(0,1), estimate.time = F, survivalOutcome = T, variance.method = "ic"), "Variance estimate is based on")
  regimes <- list(function(row) c(1,1), function(row) c(0,0))
  summary.measures <- array(1:0, dim=c(2, 1, 1))
  colnames(summary.measures) <- "A"
  r4 <- ltmleMSM(data=data, Anodes=Anodes, Lnodes=Lnodes, Cnodes=Cnodes, Ynodes=Ynodes, estimate.time = F, survivalOutcome = T, regimes=regimes, summary.measures=summary.measures, working.msm="Y~W1+W2+A", final.Ynodes = "Y2") 
  r5 <- ltmleMSM(data=data, Anodes=Anodes, Lnodes=Lnodes, Cnodes=Cnodes, Ynodes=Ynodes, estimate.time = F, survivalOutcome = T, regimes=regimes, summary.measures=summary.measures, working.msm="Y~W1+W2+A", final.Ynodes = "Y2", variance.method = "iptw") 
  expect_warning(r6 <- ltmleMSM(data=data, Anodes=Anodes, Lnodes=Lnodes, Cnodes=Cnodes, Ynodes=Ynodes, estimate.time = F, survivalOutcome = T, regimes=regimes, summary.measures=summary.measures, working.msm="Y~W1+W2+A", final.Ynodes = "Y2", variance.method = "ic"), "Variance estimate is based on") 
  
})

test_that("warning message prints when variance estimate ratio > 100", {
  n <- 1000
  W <- rnorm(n)
  A <- rexpit(8 * W)
  Y <- rexpit(W + A)
  data <- data.frame(W, A, Y)
  expect_warning(print(summary(ltmle(data=data, Anodes="A", Ynodes="Y", abar = 1, estimate.time=FALSE, gbounds = c(0, 1)))), "When this ratio is greater than 100, both variance estimates are less likely to be accurate.")
})

test_that("Robust variance is computed for non-binary outcomes", {
  n <- 500
  W1 <- rnorm(n)
  W2 <- rnorm(n)
  A1 <- rexpit(W1 + W2)
  L1 <- rnorm(n) + A1
  C1 <- rexpit(2 + W1)
  Y1 <- plogis(W1 + A1 + rnorm(n))
  A2 <- rexpit(W1 + A1)
  L2 <- rnorm(n) + L1
  C2 <- rexpit(2 + L1)
  Y2 <- 10 * plogis(W1 - L1 + A1 + A2 + rnorm(n))
  data <- data.frame(W1, W2, A1, L1, C1=BinaryToCensoring(is.uncensored=C1), Y1, A2, L2, C2=BinaryToCensoring(is.uncensored=C2), Y2)
  nodes <- list(Anodes=c("A1", "A2"), Lnodes=c("L1", "L2"), Cnodes=c("C1", "C2"), Ynodes=c("Y1", "Y2"))
  regimes <- list(function(row) c(1, 1), function(row) c(0, 0))
  summary.measures <- array(1:0, dim=c(2, 1, 1))
  colnames(summary.measures) <- "A"
  rules <- list(function(row) c(1, 1), function(row) c(as.numeric(row$W1 > 0), 1)) #regimes overlap when W1 > 0, so treatment is not static
  expect_robust_variance <- function(r) {
    expect_false(is.null(r$variance.estimate)) #NULL means variance.method was switched to "ic"
    v <- as.matrix(r$variance.estimate)
    expect_true(all(is.finite(v)))
    expect_true(all(diag(v) >= 0))
  }
  for (variance.method in c("tmle", "iptw")) {
    args <- c(list(data=data, estimate.time=FALSE, survivalOutcome=FALSE, variance.method=variance.method), nodes)
    expect_robust_variance(do.call(ltmle, c(args, list(abar=c(1, 1)))))
    expect_robust_variance(do.call(ltmle, c(args, list(rule=rules))))
    expect_robust_variance(do.call(ltmle, c(args, list(abar=c(1, 1), observation.weights=runif(n, 0.5, 2)))))
    expect_robust_variance(do.call(ltmleMSM, c(args, list(regimes=regimes, summary.measures=summary.measures, working.msm="Y ~ W1 + W2 + A", final.Ynodes="Y2"))))
  }
})

test_that("Final-node variance for non-binary outcomes is E[Y^2 | past] - Q^2, not Q(1 - Q)", {
  n <- 1000
  W <- rnorm(n)
  A <- rexpit(0.5 * W)
  Y <- plogis(W + A + rnorm(n, sd=0.1)) #in (0, 1), not binary, tiny Var(Y | W, A)
  for (variance.method in c("tmle", "iptw")) {
    r <- ltmle(data.frame(W, A, Y), Anodes="A", Ynodes="Y", abar=1, estimate.time=FALSE, variance.method=variance.method)
    #measured: about 12 with Q(1 - Q) at the final node, 0.25-0.28 with E[Y^2 | past] - Q^2
    expect_lt(r$variance.estimate / var(r$IC$tmle), 1)
  }
})

test_that("Robust variance is back-transformed for rescaled non-binary outcomes", {
  n <- 1000
  W <- rnorm(n)
  A <- rexpit(W)
  Y <- plogis(W + A + rnorm(n))
  data1 <- data.frame(W, A, Y=10 * Y + 5)
  data2 <- data.frame(W, A, Y=20 * Y + 3) #same internal [0, 1] outcome as data1
  r1 <- ltmle(data1, Anodes="A", Ynodes="Y", abar=1, estimate.time=FALSE)
  r2 <- ltmle(data2, Anodes="A", Ynodes="Y", abar=1, estimate.time=FALSE)
  r3 <- ltmle(data1, Anodes="A", Ynodes="Y", abar=1, estimate.time=FALSE, Yrange=range(data1$Y))
  expect_equal(r2$variance.estimate / r1$variance.estimate, 4, tolerance=1e-6)
  expect_equal(r3$variance.estimate, r1$variance.estimate)
  e1 <- ltmle(data1, Anodes="A", Ynodes="Y", abar=list(1, 0), estimate.time=FALSE)
  e2 <- ltmle(data2, Anodes="A", Ynodes="Y", abar=list(1, 0), estimate.time=FALSE)
  expect_false(is.null(e1$variance.estimate))
  expect_equal(summary(e2)$effect.measures$ATE$std.dev / summary(e1)$effect.measures$ATE$std.dev, 2, tolerance=1e-6)
})

test_that("Robust variance is finite and non-negative when a non-binary outcome has no noise", {
  n <- 1000
  W <- rnorm(n)
  A <- rexpit(0.5 * W)
  Y <- plogis(W + A) #Var(Y | W, A) = 0, so E[Y^2 | past] - Q^2 can be negative before bounding
  for (variance.method in c("tmle", "iptw")) {
    r <- ltmle(data.frame(W, A, Y), Anodes="A", Ynodes="Y", abar=1, estimate.time=FALSE, variance.method=variance.method)
    expect_false(is.null(r$variance.estimate))
    expect_true(is.finite(r$variance.estimate))
    expect_gte(r$variance.estimate, 0)
  }
})

test_that("TMLE based variance estimate > IC based variance estimate under positivity, non-binary outcome", {
  skip_on_cran()
  ExpectRobustStdDevUnderPositivity(function(n) runif(n))
})
