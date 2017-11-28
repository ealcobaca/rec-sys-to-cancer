bmrName = list.files(path = "../output/benchmark/", pattern = '*.RData', recursive=T)
bmrPath = paste('../output/benchmark/', bmrName, sep='')


result = c()
dataNames= c()
for (i in 1:length(bmrPath)){
#for (i in 32:34){
    obj = load(bmrPath[i])
    dataNames = c(dataNames, strsplit(strsplit(bmrPath[i], split='/')[[1]][5], split='.benchmark.RData')[[1]][1])
    result = rbind(result, c(bmr$results$data$classif.randomForest$aggr, bmr$results$data$classif.svm$aggr))
}

metadata = cbind(dataNames, data.frame(result))
colnames(metadata) = c('data.name',
    paste('classif.randomForest',names(bmr$results$data$classif.randomForest$aggr),sep='.'),
    paste('classif.svm',names(bmr$results$data$classif.randomForest$aggr),sep='.'))
rownames(metadata) = metadata[,1]
##################################################################

#Statistical testing
result.svm = c()
result.rf = c()
dataNames= c() 

for (i in 1:length(bmrPath)){
    obj = load(bmrPath[i])
    dataName = strsplit(strsplit(bmrPath[i], split='/')[[1]][5], split='.benchmark.RData')[[1]][1]
    result.rf =  bmr$results$data$classif.randomForest$measures.test$ber
    result.svm = bmr$results$data$classif.svm$measures.test$ber
    metadata[dataName,'t.test.p.value'] = t.test(result.rf, result.svm)$p.value
    metadata[dataName,'wilcox.test.p.value'] = wilcox.test(result.rf, result.svm)$p.value
}

#p-value <0.05, we can reject the null hypothesis, they are different

########################################################33333

mfName = list.files(path = "../output/meta-feature/", pattern = '*.RData', recursive=T)
mfPath = paste('../output/meta-feature/', mfName, sep='')

metadatacp = metadata
#for(i in 32:34){
for( i in 1:length(mfPath)){
    load(mfPath[i])
    dataName = strsplit(strsplit(mfPath[i], split='/')[[1]][5], split='.meta-features.RData')[[1]][1]
    obj = unlist(obj)
    for(j in names(obj)){
        metadata[as.character(dataName), j] = obj[j]
    }
}

badFeatures <- colnames(metadata)[colSums(is.na(metadata)) > 0]
ids = which(colnames(metadata) %in% badFeatures)
finalMetadata = metadata[,-ids]

write.table(file='../output/finalMetadata.dat', x=finalMetadata)

