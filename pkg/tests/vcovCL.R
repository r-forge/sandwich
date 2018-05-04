library("sandwich")
data("PetersenCL", package = "sandwich")
m <- lm(y ~ x, data = PetersenCL)

vcovCL(m, cluster = ~ firm, type = "HC0", cadjust = TRUE)
vcovCL(m, cluster = ~ firm, type = "HC0", cadjust = FALSE)
vcovCL(m, cluster = ~ firm, type = "HC1", cadjust = TRUE)
vcovCL(m, cluster = ~ firm, type = "HC1", cadjust = FALSE)
vcovCL(m, cluster = ~ firm, type = "HC2", cadjust = TRUE)
vcovCL(m, cluster = ~ firm, type = "HC2", cadjust = FALSE)
vcovCL(m, cluster = ~ firm, type = "HC3", cadjust = TRUE)
vcovCL(m, cluster = ~ firm, type = "HC3", cadjust = FALSE)

vcovCL(m, cluster = ~ firm + year, type = "HC0", cadjust = TRUE)
vcovCL(m, cluster = ~ firm + year, type = "HC0", cadjust = FALSE)
vcovCL(m, cluster = ~ firm + year, type = "HC1", cadjust = TRUE)
vcovCL(m, cluster = ~ firm + year, type = "HC1", cadjust = FALSE)
vcovCL(m, cluster = ~ firm + year, type = "HC2", cadjust = TRUE)
vcovCL(m, cluster = ~ firm + year, type = "HC2", cadjust = FALSE)
vcovCL(m, cluster = ~ firm + year, type = "HC3", cadjust = TRUE)
vcovCL(m, cluster = ~ firm + year, type = "HC3", cadjust = FALSE)


# vcovCL compared with cluster.vcov from multiwayvcov package

(cl1 <- vcovCL(m, cluster = ~ firm, type = "HC1", cadjust = TRUE))
(cl2 <- multiwayvcov::cluster.vcov(m, ~ firm))

(cl3 <- vcovCL(m, cluster = ~ firm + year, multi0 = TRUE))
(cl4 <- multiwayvcov::cluster.vcov(m, cbind(PetersenCL$firm, PetersenCL$year)))

(cl5 <- vcovCL(m, cluster = ~ firm, type = "HC0", cadjust = FALSE))
(cl6 <- multiwayvcov::cluster.vcov(m, ~ firm, df_correction = FALSE))

(cl7 <- vcovCL(m, cluster = ~ firm + year, type = "HC0", cadjust = FALSE))
(cl8 <- multiwayvcov::cluster.vcov(m, cbind(PetersenCL$firm, PetersenCL$year), df_correction = FALSE))

all.equal(cl1, cl2)
all.equal(cl3, cl4)
all.equal(cl5, cl6)
all.equal(cl7, cl8)


# vcovCL compared with BMlmSE (Bell-McCaffrey standard errors as described in Imbens and Kolesar 2016,
# https://github.com/kolesarm/Robust-Small-Sample-Standard-Errors)
# BMlmSE(m, clustervar = factor(PetersenCL$firm), IK = FALSE)

bellmc1 <- matrix(c(4.494487e-03, -6.592912e-05, -6.592912e-05, 2.568236e-03), nrow = 2)
rownames(bellmc1) <- colnames(bellmc1) <- c("(Intercept)", "x")
bellmc1
(bellmc2 <- vcovCL(m, cluster = ~ firm, type = "HC2", cadjust = TRUE))

all.equal(bellmc1, bellmc2, tolerance = 1e-7)


data("InstInnovation", package = "sandwich")
n <- glm(cites ~ institutions, family = "poisson", data = InstInnovation)

vcovCL(n, cluster = ~ company, type = "HC0", cadjust = TRUE)
vcovCL(n, cluster = ~ company, type = "HC0", cadjust = FALSE)
vcovCL(n, cluster = ~ company, type = "HC1", cadjust = TRUE)
vcovCL(n, cluster = ~ company, type = "HC1", cadjust = FALSE)
vcovCL(n, cluster = ~ company, type = "HC2", cadjust = TRUE)
vcovCL(n, cluster = ~ company, type = "HC2", cadjust = FALSE)
vcovCL(n, cluster = ~ company, type = "HC3", cadjust = TRUE)
vcovCL(n, cluster = ~ company, type = "HC3", cadjust = FALSE)

vcovCL(n, cluster = ~ company + year, type = "HC0", cadjust = TRUE)
vcovCL(n, cluster = ~ company + year, type = "HC0", cadjust = FALSE)
vcovCL(n, cluster = ~ company + year, type = "HC1", cadjust = TRUE)
vcovCL(n, cluster = ~ company + year, type = "HC1", cadjust = FALSE)
vcovCL(n, cluster = ~ company + year, type = "HC2", cadjust = TRUE)
vcovCL(n, cluster = ~ company + year, type = "HC2", cadjust = FALSE)
vcovCL(n, cluster = ~ company + year, type = "HC3", cadjust = TRUE)
vcovCL(n, cluster = ~ company + year, type = "HC3", cadjust = FALSE)


o <- lm(log(cites) ~ institutions, data = InstInnovation, subset = cites > 0)

vcovCL(o, cluster = ~ company + year, type = "HC0", cadjust = TRUE, multi0 = TRUE)
vcovCL(o, cluster = ~ company + year, type = "HC0", cadjust = TRUE, multi0 = FALSE)
vcovCL(o, cluster = ~ company + year, type = "HC0", cadjust = FALSE, multi0 = TRUE)
vcovCL(o, cluster = ~ company + year, type = "HC0", cadjust = FALSE, multi0 = FALSE)
vcovCL(o, cluster = ~ company + year, type = "HC1", cadjust = TRUE, multi0 = TRUE)
vcovCL(o, cluster = ~ company + year, type = "HC1", cadjust = TRUE, multi0 = FALSE)
vcovCL(o, cluster = ~ company + year, type = "HC1", cadjust = FALSE, multi0 = TRUE)
vcovCL(o, cluster = ~ company + year, type = "HC1", cadjust = FALSE, multi0 = FALSE)
vcovCL(o, cluster = ~ company + year, type = "HC2", cadjust = TRUE, multi0 = TRUE)
vcovCL(o, cluster = ~ company + year, type = "HC2", cadjust = TRUE, multi0 = FALSE)
vcovCL(o, cluster = ~ company + year, type = "HC2", cadjust = FALSE, multi0 = TRUE)
vcovCL(o, cluster = ~ company + year, type = "HC2", cadjust = FALSE, multi0 = FALSE)
vcovCL(o, cluster = ~ company + year, type = "HC3", cadjust = TRUE, multi0 = TRUE)
vcovCL(o, cluster = ~ company + year, type = "HC3", cadjust = TRUE, multi0 = FALSE)
vcovCL(o, cluster = ~ company + year, type = "HC3", cadjust = FALSE, multi0 = TRUE)
vcovCL(o, cluster = ~ company + year, type = "HC3", cadjust = FALSE, multi0 = FALSE)



