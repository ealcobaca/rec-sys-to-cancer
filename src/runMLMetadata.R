library('FSelector')
library(mlr)
library("ggthemes")
library("scales")
library(ggplot2)

seed=123
set.seed(seed = seed)
options(mlr.debug.seed = seed)


bmr = function(task){
    lrns = list( 
                classif.randomForest=makeLearner("classif.randomForest"),
                classif.svm=makeLearner("classif.svm", kernel = "linear"),
                classif.C50 = makeLearner("classif.C50"),
                classif.knn = makeLearner("classif.knn"),
                classif.naiveBayes = makeLearner("classif.naiveBayes")
                )

    rdesc = makeResampleDesc("LOO")
    bmr = benchmark(lrns, task, rdesc, 
                    measures=list(acc, mmce, ber, kappa, timetrain), show.info=F)

    return(bmr)

}

all.data = function(){
    md = read.table('../output/finalMetadata.dat')
    ### rule to classification
    Class = ifelse( md$classif.randomForest.ber.test.mean < md$classif.svm.ber.test.mean, 
                   'classif.randomForest', 'classif.svm')
    metadata = Filter(function(x) sd(x) != 0, md[,-c(1:11)])
    metadata = cbind(metadata, Class)

    types = c('^general','^statistical','^infotheo','^model.based', '^landmarking')
    classId = ncol(metadata)
    result = list()
    count = 1

    print('all')
    task = makeClassifTask(data = metadata, target = "Class")
    result[[count]]= bmr(task)

    for (type in types){
        print(type)
        count = count +1
        ids = c(grep(type,colnames(metadata)), classId)

        ### models generator
        task = makeClassifTask(data = metadata[,ids], target = "Class")
        result[[count]]= bmr(task)


    }


    types = c('pca.10','pca.25','pca.40','pca.55')
    for (i in c(10,25,40,55)){
        aux1 = metadata
        metadata = cbind(prcomp(metadata[,-ncol(metadata)], center = TRUE, scale. = TRUE, rank.=i)$x, 
            data.frame(Class=as.character(metadata$Class)))

        count = count +1

        ### models generator
        task = makeClassifTask(data = metadata, target = "Class")
        result[[count]]= bmr(task)
        metadata = aux1

    }

    weights <- relief(Class~., metadata, neighbours.count = 5, sample.size = 20)
    #print(weights)

    for (i in c(10,25,40,55)){
        count = count + 1
        subset <- c(cutoff.k(weights, i),'Class')
        df = metadata[,subset]
    
        task = makeClassifTask(data = metadata, target = "Class")
        result[[count]]= bmr(task)
    }

    names(result) = c('all', 'general','statistical','infotheo','model.based', 'landmarking',
                      'pca.10','pca.25','pca.40','pca.55', 'relief.10', 'relief.25', 'relief.40', 'relief.55')

    return(result) 

}

agg.bmrs = function(bmrs){
    
    md = read.table('../output/finalMetadata.dat')
    ### rule to classification
    Class = ifelse( md$classif.randomForest.ber.test.mean < md$classif.svm.ber.test.mean,
        'classif.randomForest', 'classif.svm')

    baseline.random = mean(sample(c(0,1), 86, replace=T))
    baseline.majority.class = 0.5
    baseline.tech.diff = 29/length(Class)


    tp = names(bmrs)
    i = 0
    df.perf.data = data.frame(perf=NA, classif=NA, tp.mf=NA)
    for ( bmr in bmrs){
        i = i+1
        perf = c(bmr$results$metadata$classif.randomForest$aggr[1],
                 bmr$results$metadata$classif.svm$aggr[1],
                 bmr$results$metadata$classif.C50$aggr[1],
                 bmr$results$metadata$classif.knn$aggr[1],
                 bmr$results$metadata$classif.naiveBayes$aggr[1],
#                 bmr$results$metadata$classif.xgboost$aggr[1],
                 baseline.random,
                 baseline.majority.class
                 )

        classif = c('classif.randomForest',
                    'classif.svm',
                    'classif.C50',
                    'classif.knn',
                    'classif.naiveBayes',
#                    'classif.xgboost',
                    'baseline.random',
                    'baseline.majority.class')

        tp.mf = rep(tp[i],7)
        aux = data.frame(perf=perf,classif=classif, tp.mf=tp.mf)
        df.perf.data = rbind(df.perf.data, aux)

    }


    return(df.perf.data[-1,])
}

heat.map = function(df.perf.data, title, path){
    
    l = unique(df.perf.data$tp.mf)
    m = c()
    for(i in l){
        m = c(m, mean(df.perf.data[df.perf.data$tp.mf == i, ]$perf))
    }
    l = l[order(m)]
    df.perf.data$tp.mf = factor(df.perf.data$tp.mf, levels=l)

    l = unique(df.perf.data$classif)
    m = c()
    for(i in l){
        m = c(m, mean(df.perf.data[df.perf.data$classif == i, ]$perf))
    }
    l = l[order(m)]
    l = c("baseline.majority.class", "baseline.random",l[c(-which(l=="baseline.majority.class"), -which(l=="baseline.random"))])
    df.perf.data$classif = factor(df.perf.data$classif, levels=l)

    myPlot=    ggplot(data = df.perf.data, aes(x = classif, y = tp.mf)) +
        geom_tile(aes(fill = perf))+
        scale_fill_gradient2(name = "Acur\u{00E1}cia", low = "blue", mid = "white", high = "Red",
            limits = c(0.43, 0.85), midpoint = 0.65, breaks=c(0.45,0.55,0.65,0.75,0.85))+
        geom_text(aes(label = round(perf, 2)), fontface='bold') +
        labs(title=title, x='Classificadores', y='Meta-Bases') +
        theme(plot.title = element_text(hjust = 0.5, face='bold'), axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"), legend.text=element_text(size=14), legend.title=element_text(size=14, hjust = -1),
        axis.text.x = element_text(angle = 90, hjust = 1))


    ggsave(filename=path, plot=myPlot, width=14, height=7, dpi=300)

}


rf.eval = function(result){
   
    model = result$general$results$'metadata[, ids]'$classif.randomForest$model[[1]]
    imp = model$learner.model$importance
    imp = data.frame(imp)
    imp$Features = rownames(imp)

    imp = imp[order(imp$MeanDecreaseGini,decreasing=TRUE),]
    imp$Features = factor(imp$Features, levels=unique(imp$Features))

    p = ggplot(data=imp, aes(x=Features, y=MeanDecreaseGini)) +
        geom_bar(stat="identity")+        theme(plot.title = element_text(hjust = 0.5, face='bold'), axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1))+
        ggtitle("Import\u00e2ncia do RF")



    ggsave(filename='../output/imgs/imp.pdf', plot=p, width=6, height=6, dpi=300)


    error = result$general$results$'metadata[, ids]'$classif.randomForest$pred$data[,2:3]
    error = ifelse(error[,1] != error[,2], 'Erro', 'Acerto')

     metadata = read.table('../output/finalMetadata.dat')
    data = metadata[,1:11]
    data$data.name = c(rep('RNA-Seq',length(grep('rnaseqv2',data$data.name))),
                       rep('miRNA-Seq',length(grep('mirnaseq',data$data.name))),
                       rep('microarray',length(grep('microarray',data$data.name))))
    aux = data[,c('data.name','classif.svm.ber.test.mean', 'classif.randomForest.ber.test.mean', 'wilcox.test.p.value')]
    aux$diff = (1-aux$classif.svm.ber.test.mean) - (1-aux$classif.randomForest.ber.test.mean)
    aux = aux[order(aux[,c('diff')]),]
    aux$id = 1:86
    aux$error = factor(error, levels=c('Erro','Acerto'))

     p1=ggplot(aux, aes(x=id, y=diff, fill=data.name,alpha=error)) +
        xlab("Conjunto de Dados") +
        ylab("Diferen\u{00E7}a BAC (SVM - RF)") +
        ggtitle("Diferen\u{00E7}a de desempenho (BAC) entre SVM e RF") +
        geom_bar(stat = 'identity')+
        scale_color_tableau(palette="colorblind10")+
        scale_fill_tableau(palette="colorblind10")+
        scale_alpha_discrete(range=c(0.4,1))+
        annotate("text", x = 75, y =0.35, label = "Ganho com SVM", size=5 , fontface="bold")+
        annotate("text", x = 15, y =-0.35, label = "Ganho com RF", size=5 , fontface="bold")+
        guides(fill=guide_legend(title="Tipo do Conjunto de Dado"), alpha=guide_legend(title='Transpar\u{00EA}ncia'))+
        theme(plot.title = element_text(hjust = 0.5, face='bold'), axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"))



    multiplot(p1, cols=1)

}

result = all.data()
df.perf.data = agg.bmrs(result)
heat.map(df.perf.data[1:42,], 'Heat Map do LOOCV com diferentes conjuntos de meta-caracteristicas','../output/imgs/diff.pdf')
heat.map(df.perf.data[43:70,], "Heat Map do LOOCV aplicando PCA em todas as meta-caracteristicas",'../output/imgs/pca.pdf')
heat.map(df.perf.data[71:98,], "Heat Map do LOOCV aplicando RELIEF em todas as meta-caracteristicas",'../output/imgs/fsel.pdf')

