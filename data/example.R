#
# Prediction of relatedness between genes under given certain condictions with 
# MFR model

#
# Dependence:
# R package ¡®e1071¡¯: an R language version of LIBSVM.

# 
# Input:

# x: a data.frame object, with 12 columns, each of which is the similarity-
# based gene pair features, including expression level of the first gene(exp1), 
# expression level of the second gene (exp2), co-expression levels calculatded 
# using Pearson correlation coefficient (PCC), Spearman rank correlation (SRC), 
# mutual information (MI), partial Pearson correlation (PPC), conditional mutual 
# information (CMI), Gene Ontology similarity (goSim), Reactome similarity 
# (rxSim), subcellular localization similarity (lcSmi), homology similarity 
# (hgSim), and transcriptional regulatory similarity (trSim). For details, see 
# dataDemo.Rdata.
 
# model: a trained SVM based MFR model.

#
# Output:
  
# MFR values of gene pairs. The larger value indicates the two genes are more 
# considerably interacting.

#
# function:

predict.mfr = function(x, model = model.mfr){
    require(e1071)
    data.dn = x
    sv = model$SV
    cs = model$coefs
    w = t(sv) %*% cs
    d = model$rho
    x1= matrix(0, nrow(data.dn), ncol(data.dn))
    for(i in c(1:12)){
        x1[,i] = as.numeric(data.dn[,i]) 
    }
    y1 = matrix(model$x.scale$`scaled:center`, 
                nrow(data.dn), ncol(data.dn), byrow = T)
    s1 = matrix(model$x.scale$`scaled:scale`, 
                nrow(data.dn), ncol(data.dn), byrow = T)
    xx = (x1 - y1) / s1
    
    mfr = xx %*% w - d
    mfr = mfr * model$y.scale$`scaled:scale` + model$y.scale$`scaled:center`
    return(sigmoid(mfr))
}

#
#an example:

# ... is the path of the resources
load(".../model.mfr.RData")
load(".../data.demo.RData")

mfr.values = predict.mfr(x = data.x)
print(mfr.values)





