library("sandwich")
data("PetersenCL", package = "sandwich")
m <- lm(y ~ x, data = PetersenCL)

vcovPL(m, cluster = PetersenCL$firm, adjust = TRUE)
vcovPL(m, cluster = PetersenCL$firm, adjust = FALSE)

data("InstInnovation", package = "sandwich")
n <- glm(cites ~ institutions, family = "poisson", data = InstInnovation)

vcovPL(n, cluster = InstInnovation$industry, adjust = TRUE)
vcovPL(n, cluster = InstInnovation$industry, adjust = FALSE)


