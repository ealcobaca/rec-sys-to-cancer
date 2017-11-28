library(mlr)

set.seed(123)

#task = makeClassifTask(data = iris, target = "Species")
#classif.rf = makeLearner("classif.randomForest")
#classif.svmL = makeLearner("classif.svm", kernel = "linear")
#rdesc = makeResampleDesc("RepCV", reps = 10)

## Calculate the performance measures
#result.rf = resample(classif.rf, task, rdesc, measures = list(mmce, timetrain))
#result.svm = resample(classif.rf, task, rdesc, measures = list(mmce, timetrain))

dataName = list.files(path = "../datasets/pca/", pattern = '*.csv', recursive=T)
dataPath = paste('../datasets/pca/', dataName, sep='')

classif.rf = makeLearner("classif.randomForest")
classif.svmL = makeLearner("classif.svm", kernel = "linear")

## Two learners to be compared
lrns = list(classif.rf,classif.svmL)



result = list()
i=0
for(f in dataPath){
    i = i+1
    cat('interation [',i,']\n')
    cat('file = ',f,'\n\n')
    data = read.csv2(f)
    data$Class = factor(data$Class) 
    #print(data$Class)
    task = makeClassifTask(data = data, target = "Class")
    
    ## Choose the resampling strategy
    rdesc = makeResampleDesc("RepCV", reps = 10)

    ## Conduct the benchmark experiment
    bmr = benchmark(lrns, task, rdesc, 
        measures=list(acc, mmce, ber, kappa, timetrain), show.info=FALSE)

    save(bmr,file=paste('../output/benchmark/',
        paste(strsplit(dataName[i],'.csv')[[1]],
              '.benchmark.RData',sep='')
              ,sep=''))

}

############################################

