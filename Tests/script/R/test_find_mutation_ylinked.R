source("mutation_function.R")
options(digits=10)


## test_prior
theta<- 0.001
freq<- c(0.3, 0.2, 0.2, 0.3)
ref_weight<- 1

freq_sum<- sum(freq)
normalised_theta<- theta/freq_sum

genotype_prior<- vector(length=5,mode="list")

for(i in 1:5){
    alpha<- freq * normalised_theta
    if(i<5){
        alpha[i] <- alpha[i] + ref_weight
    }
    alpha_sum<- sum(alpha)
    genotype_prior[[i]] <- alpha/alpha_sum
}    



for(i in 1:5){
    cat("{", paste(genotype_prior[[i]], collapse=", "), "},", sep="", fill=T)
}

