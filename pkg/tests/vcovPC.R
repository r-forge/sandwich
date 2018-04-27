library("sandwich")
data("PetersenCL", package = "sandwich")
m <- lm(y ~ x, data = PetersenCL)

vcovPC(m, cluster = PetersenCL$firm, order.by = PetersenCL$year)

PU <- subset(PetersenCL, !(firm == 1 & year == 10))
u_m <- lm(y ~ x, data = PU)

vcovPC(u_m, cluster = PU$firm, order.by = PU$year, pairwise = TRUE)
vcovPC(u_m, cluster = PU$firm, order.by = PU$year, pairwise = FALSE)


