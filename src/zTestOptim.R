source("src/mutation.R")
options(digits=10)

rd<-vector(length=3, mode="list")
genotype<-vector(length=3, mode="list")
rd[[1]]<-c(0, 1, 25, 29) #878
genotype[[1]]<-c(1.9681812586082493465e-19, 5.8848619632386626597e-17, 9.3761538576728770258e-11, 9.3067830049776481198e-11, 1.1750042113891158116e-16, 2.571990828847844964e-08, 2.7827281184883153803e-08, 1.817757848706868286e-07, 1, 8.3192969181890452131e-07)

rd[[2]]<-c(0, 0, 57, 0) #891
genotype[[2]]<-c(2.9503055407498715637e-18, 2.9503055407498715637e-18, 2.2779180406756631894e-07, 2.9503055407498715637e-18, 2.9503055407498715637e-18, 2.2779180406756631894e-07, 2.9503055407498715637e-18, 1, 2.2779180406756631894e-07, 2.9503055407498715637e-18)

rd[[3]]<-c(0, 0, 76, 1) #892
genotype[[3]]<-c(1.6296875623801300125e-19, 1.6296875623801300125e-19, 7.848024187308119499e-08, 4.8727658115165870234e-17, 1.6296875623801300125e-19, 7.848024187308119499e-08, 4.8727658115165870234e-17, 1, 2.152806634867116697e-05, 9.7292347474093002355e-17)

upper<-vector(length=2, mode="list")  #1 dad, 2 mom
upper[[1]]<- c(0.00014982019479770606772, 5.9910104887616127319e-08, 0.00029961043454296829892, 8.9865157331424197595e-08, 9.9870144847656117711e-05, 0.00019974028969531223542, 5.9910104887616127319e-08, 0.99880131862140886234, 0.00029961043454296829892, 0.00014982019479770606772)
upper[[2]]<- upper[[1]]

lower<-vector(length=2, mode="list") #1 dad, 2 mom

freq<-c(0.3, 0.2, 0.2, 0.3)
F81_mu3<- f81_full(3e-8,freq)
F81_mu2<- f81_full(2e-8,freq)

numMutation<- 0; #0, 1, -1, 2
mut_somatic<- mitosis_diploid_matrix(F81_mu2, numMutation)

lower[[1]]<- mut_somatic %*% genotype[[1]]
lower[[2]]<- mut_somatic %*% genotype[[2]]
resultNo<- sum(lower[[2]] * lower[[1]])

numMutation<- 0; #0, 1, -1, 2
mut_somatic<- mitosis_diploid_matrix(F81_mu2, -1)

lower[[1]]<- mut_somatic %*% genotype[[1]]
lower[[2]]<- mut_somatic %*% genotype[[2]]
resultFull<- sum(lower[[2]] * lower[[1]])

resultNo
resultFull

calLikeli <- function(branch, genotype){
  rate<- branch#*1e-6
  if(any(rate > 1 |  rate < 1e-30 )){
    return(-Inf)
  }
  # if(any(branch > 100000 |  branch < 1e-50 )){
  #   return(-Inf)
  # }

  mut1<- mitosis_diploid_matrix(f81_full(rate[1], freq), -1)
  mut2<- mitosis_diploid_matrix(f81_full(rate[2], freq), -1)

  r<- (mut1 %*% genotype[[1]]) * (mut2 %*% genotype[[2]])
  resultFull<- sum(log(apply(r,2,sum)))
  return(resultFull)

}
# optimize(calLikeli, c(0.1,1e-100), genotype=genotype, maximum=T)



genotype<- list()
nSite<- 1e6
genotype[[1]]<- matrix(runif(nSite*10), nrow=10)
genotype[[2]]<- matrix(runif(nSite*10), nrow=10)
# genotype[[1]]<- matrix(rep(0, nSite*10), nrow=10)
# genotype[[2]]<- matrix(rep(0, nSite*10), nrow=10)



G0<- seq(0,(nSite-1)*10,by=10)+  floor(runif(nSite,1,11))

key<- floor(runif(nSite,1,11))
G1<- seq(0,(nSite-1)*10,by=10)+ key
G2<- G1
gmap<- c(2,3,4,2,2,3,5,6,8,9)
x<- sample(nSite,nSite*0.1025)
G2[x]<- seq(0,(nSite-1)*10,by=10)[x] + gmap[ key[x] ]

# genotype[[1]][G0]<- 1
# genotype[[2]][G0]<- 1
genotype[[1]][G1]<- genotype[[1]][G1]+100
genotype[[2]][G2]<- genotype[[2]][G2]+100


# genotype[[2]]<- genotype[[1]] #+runif(10)
genotype<- lapply(genotype, function(x){
  apply(x, 2, function(y){
    y/sum(y)
  })
})


# branch<- c(0.1,10)
# optim(rep(runif(1),2), calLikeli, genotype=genotype[1:2], control=list(fnscale=-1))$par
branch<- runif(2,0,1)
optim(branch, calLikeli, genotype=genotype, control=list(fnscale=-1
  # ,reltol=1e-16
 # ,trace=T
 ,maxit=1000
));branch


branch<- c(0.664391039405018, 0.379738186253235)

$par
[1] 0.02413553142192576 0.00536953347295528

$par
[1] 24152.1284130263  5347.5723149858

$value
[1] -530.004221502977
calLikeli(branch, genotype)

optim(c(1e-3,1e-8), calLikeli, genotype=genotype[1:2], control=list(fnscale=-1))
optim(c(1e-8,1e-3), calLikeli, genotype=genotype[1:2], control=list(fnscale=-1))
optim(c(0.019,0.099), calLikeli, genotype=genotype[1:2], control=list(fnscale=-1))




a=apply(genotype[[1]],2,which.max)
b=apply(genotype[[2]],2,which.max)
plot(a-b)
plot(abs(a-b))


mut1<- matrix(c(0.9,0.2,0.1,0.8),nrow=2, byrow=T)
mut2<- matrix(c(0.6,0.4,0.3,0.7),nrow=2, byrow=T)
gg<- c(0.9,0.1)

sum( mut1 %*% gg * mut1 %*% gg )
sum( mut2 %*% gg * mut2 %*% gg )
