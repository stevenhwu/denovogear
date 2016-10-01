source("mutation_function.R")
options(digits=10)


## test_prior
theta<- 0.001
freq<- c(0.3, 0.2, 0.2, 0.3)
ref_weight<- 1

sum_freq<- sum(freq)
normalised_theta<- theta/sum_freq

genotype_prior<- vector(length=5,mode="list")
for(i in 1:5){


    genotype_prior[[i]]<- freq * normalised_theta
    if(i<5){
        genotype_prior[[i]][i] <- genotype_prior[[i]][i] + ref_weight
    }
    



}
