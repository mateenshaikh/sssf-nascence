rm(list=ls())
gc()

set.seed(123)
n=10
trueb = c(.25)
x = runif(n)

residuals = rnorm(n,0,sd=.5)
y = x*trueb + residuals

ls3 = function(betavals,k){
  sapply(betavals,function(beta)
  n*log(sum((y-beta*x)^2)) +1 -1/(cosh(beta*k))
  )
}

b = seq(-.05,.35,.001)
kvals = c(100,20,10,5,1)
lsvalues = sapply(kvals,function(k)sapply(b,ls3,k=k))

g = ((1:length(kvals))/(length(kvals)*2))
roots = numeric(0)
xx = sum(x^2)
intervals = list()
for (i in 1:length(kvals)){
  k = kvals[i]
  f = function(bb)  sapply(bb,function(beta){
    res = y - x*beta
    res2 = sum(res^2)
    -2*n/res2*sum(res*x) + k*(tanh(beta*k))/(cosh(beta*k))
  })
  initialspot = 0
  if(i>1)initialspot = roots[i-1]
  interv = initialspot + c(-.01,.01)
  intervals[[i]] = interv
  # curve(f,.7,interv[2])
  roots[i] = uniroot(f,interval=interv,extendInt = "yes",tol = 1e-6)$root
} 

fvalues = sapply(1:length(roots),function(i)ls3(b=roots[i],k=kvals[i]))

suppressPackageStartupMessages(require(tidyverse))
opt_df = data.frame(
  i = seq_along(kvals),
  k = kvals,
  beta_opt = roots
) |>
  mutate(aic_opt = map2_dbl(beta_opt, k, ls3))

# ---- curves ----
curves_df = expand.grid(beta = b, k = kvals) |>  mutate(aic = map2_dbl(beta, k, ls3))

ref_k = 10^5
ref_df = data.frame(beta = b) |>
  mutate(aic = map_dbl(beta, \(x) ls3(x, ref_k)))

i_focus = 1
k_focus = kvals[i_focus]

curves_focus_df = curves_df |> filter(k == k_focus)

opt_focus_df = opt_df |>
  filter(i >= i_focus) |>
  mutate(
    label = paste0("k=", k),
    dx = .75*if_else(row_number() == 1, -0.03, 0.01),
    dy = .75*if_else(row_number() == 1, -0.10, 0.10)
  )

# ---- plot ----
pdf(file="aicplotexample.pdf",width = 8,height=8)
ggplot() +
  geom_line(
    data = curves_df,
    aes(x = beta, y = aic, group = k),
    color = "grey30"
  ) +
  geom_line(
    data = ref_df,
    aes(x = beta, y = aic),
    color = "red",
    linewidth = 0.7
  ) +
  geom_line(
    data = opt_focus_df,
    aes(x = beta_opt, y = aic_opt),
    color = "black",
    linewidth = 0.6,
    linetype = 2
  ) +
  geom_point(
    data = opt_focus_df,
    aes(x = beta_opt, y = aic_opt),
    color = "black",
    size = 2
  ) +
  geom_text(
    data = opt_focus_df,
    aes(x = beta_opt + dx, y = aic_opt + dy, label = label),
    size = 6,
    color = "black"
  ) +
  labs(
    x = expression(beta),
    y = "Approximated AIC",
    title = "Sequence of surrogate functions and their optima"
  ) +
  coord_cartesian(ylim = range(curves_df$aic)) +
  theme_bw()+
  theme(text = element_text(size = 20))


dev.off()
