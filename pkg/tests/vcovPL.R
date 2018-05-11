library("sandwich")
data("PetersenCL", package = "sandwich")
m <- lm(y ~ x, data = PetersenCL)

vcovPL(m, cluster = PetersenCL$firm, adjust = TRUE)
vcovPL(m, cluster = PetersenCL$firm, adjust = FALSE)

data("InstInnovation", package = "sandwich")
n <- glm(cites ~ institutions, family = "poisson", data = InstInnovation)

vcovPL(n, cluster = InstInnovation$industry, adjust = TRUE)
vcovPL(n, cluster = InstInnovation$industry, adjust = FALSE)


# vcovPL covariance matrix compared with vcovSCC from plm package
pm <- plm::plm(y ~ x, data = PetersenCL, model = "pooling", indexes = c("firm", "year"))

(pl1 <- vcovPL(m, cluster = PetersenCL$firm, adjust = FALSE))
(pl2 <- plm::vcovSCC(pm))
attr(pl2, "cluster") <- NULL

all.equal(pl1, pl2)


# vcovPL compared with Stata's xtscc (xtscc y x, lag(1) ase) standard errors
pl4 <- c(0.0243573, 0.0281633)
names(pl4) <- c("(Intercept)", "x")
pl4
(pl3 <- sqrt(diag(pl1)))

all.equal(pl3, pl4, tolerance = 1e-6)
