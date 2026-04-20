setwd("~/Dropbox/tru/research/spl/noclass/rsim")
rm(list = ls())
bmatgaussian_r = function(x.orig,logkseq,numiter = 1,SAFETY = 0,onlyfinalsols = F,icpenalty = 2){
  z = order(x.orig)
  x = x.orig[z] #note that x = x[order(order(x))]
  n = length(x)
  ic = icpenalty/2
  muhat = x
  sigma2=1
  bmat = matrix(NA,nrow = n,ncol = length(logkseq))
  
  bic1 = bic2 = bic3 = 0
  for (kindex in 1:length(logkseq)){
    k = 2^logkseq[kindex]
    k2 = k^2
    
    for (i in 1:n){
      distances = (muhat-muhat[i])
      weights = 1/cosh(k*distances)#weight = "distance on 0-1 scale"

      bic1 = sum((x-muhat)/(sigma2));
      bic2 = -n;
      bic3 = 0;
      
      
      if (i>1){
        th = tanh(k*(muhat[i]-muhat[i-1]));
        sh = 1/cosh(k*(muhat[i]-muhat[i-1]));
        bic1 = bic1-ic*k*th*sh;
        bic2 = bic2-ic*k2*sh*(sh*sh-th*th);
        #          bic3 = bic3+logn*k3*sh*th*(th*th-5*sh*sh);
      }
      if (i<n){
        th = tanh(k*(muhat[i+1]-muhat[i]));
        sh = 1/cosh(k*(muhat[i+1]-muhat[i]));
        bic1 = bic1+ic*k*th*sh;
        bic2 = bic2-ic*k2*sh*(sh*sh-th*th);
        #         bic3 = bic3-logn*k3*sh*th*(th*th-5*sh*sh);
      }
      bic12 = bic1  /   (sign(bic2)*(SAFETY+abs(bic2)));
      muhat[i] = muhat[i]-bic12;
      #        muhat[i] = muhat[i]-bic3/2*(bic12)^2;
      #        muhat[i] = truemeans[i]
      
    }
    bmat[,kindex] = muhat
  }
  
  objfunval = -2*sum(dnorm(x,muhat,sqrt(sigma2),log = T))+ic*2*sum(weights)
  
  
  if (onlyfinalsols){ 
    return (list(finalsols = (bmat[,ncol(bmat)]),optval = objfunval))
  }

  bmat = bmat[order(z),  ]
  
  list (solpath = t(bmat),optval = objfunval)
}

# give each group approx unit variance by 
# giving unit variance then multiply by two
# then subset to every fifth observation
indices = seq(1,nrow(faithful),5)
scaleAndSubset = \(x) (2*scale(x))[indices]
faithful_modified = apply(faithful,2,scaleAndSubset) |> 
  as.data.frame()

logkseq = seq(-3,7,.01)
ans = apply(faithful_modified,2,\(x)bmatgaussian_r(x,logkseq),simplify = F)


makegroups = function(u,v){
  makegroup = function(a){
    y = a$solpath[nrow(a$solpath),] 
    nearmean1 = (abs(y-y[1])<1)
    nearmean1+1
  }
  y1 = makegroup(u)
  y2 = makegroup(v)
  combined = y1
  combined[y1>y2] = 4
  combined[y1<y2] = 3
  combined
}
groupings = with(ans,makegroups(eruptions,waiting))

par(mfrow = c(1,3))
matplot(logkseq,ans$eruptions$solpath,type = "l",lty = 1);
matplot(logkseq,ans$waiting$solpath,type = "l",lty = 1);
plot(faithful_modified,col = groupings,pch = groupings)


library(tidyverse)

kseq = 2^logkseq

df_ans = function(ans){
  as.data.frame(ans$solpath) |>
  mutate(kseq = kseq) |>
  pivot_longer(cols = -kseq, names_to = "variable", values_to = "value")
}

df_faithful = lapply(ans,df_ans)

# Plot
eruptions = ggplot(df_faithful$eruptions, aes(x = kseq, y = value, group = variable)) +
  geom_line(linetype = 1,linewidth = .1) + 
  theme_minimal() +scale_color_brewer(palette = "Dark2") +
  labs(title = "Solution Path for Eruption Time", x = "k", y = expression(paste(mu[eruptions])))+ggplot2::scale_x_log10()

waiting = ggplot(df_faithful$waiting, aes(x = kseq, y = value, group = variable)) +
  geom_line(linetype = 1,linewidth = .1) + 
  theme_minimal() +scale_color_brewer(palette = "Dark2") +
  labs(title = "Solution Path for Waiting Time", x = "k", y = expression(paste(mu[waiting])))+ggplot2::scale_x_log10()


groupings[groupings = 3] = "Split"

df_faithful = as.data.frame(faithful[indices, ]) |>
  mutate(Group = as.factor(groupings)) # Ensure grouping is a factor for discrete scales

cb_palette = c("#E69F00", "#0072B2",  "#009E73", "#F0E442","#56B4E9",  "#D55E00", "#CC79A7")
clustering = ggplot(df_faithful, aes(x = eruptions, y = waiting, color = Group, shape = Group)) +
  geom_point() +
  theme_minimal() +scale_color_manual(values = cb_palette) +
  labs(title = "Old Faithful Geyser Clustering", x = "Eruptions", y = "Waiting Time")


h = 8
w = 8
pdf(file = "eruptions.pdf",width = w,height = h)
eruptions
dev.off()

pdf(file = "waiting.pdf",width = w,height = h)
waiting
dev.off()

pdf(file = "clustering.pdf",width = w*1.33,height = h)
clustering
dev.off()

pdf(file = "all3.pdf",width = 8,height = 2)
library(patchwork)
eruptions+waiting+clustering
dev.off()


##comparison not implemented
# library(mclust)
# result = Mclust(faithful[indices,],G = 2:3,modelNames = "EII")
# colour = round(result$z)+1
# plot(faithful[indices,],col = colour)
