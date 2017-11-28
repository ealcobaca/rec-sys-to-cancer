library("ggplot2")
library("ggthemes")
library("scales")

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, cols=1, layout=NULL, save=FALSE, file='data.pdf') {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }

    if (save == TRUE){
        ggsave(file)    
    }
}

pca.2comp.diff = function(){
    df.pca.comp.microarray = read.csv2('../output/percentage100PcaComponents.microarray.csv')
    df.pca.comp.mirnaseq = read.csv2('../output/percentage100PcaComponents.mirnaseq.csv')
    df.pca.comp.rnaseqv2 = read.csv2('../output/percentage100PcaComponents.rnaseqv2.csv')
    mic = cbind( name=as.character(df.pca.comp.microarray[,1]), data.frame(sum=rowSums(df.pca.comp.microarray[,c(2,3)])))
    mic = mic[order(mic[,2]),]
    porWorst = df.pca.comp.microarray[df.pca.comp.microarray[,1] == as.character(mic[1,1]),2:3] * 100
    porBest = df.pca.comp.microarray[ df.pca.comp.microarray[,1] == as.character(mic[nrow(mic),1]),2:3] *100

    best= read.csv2(paste('../datasets/pca/microarray/',mic[nrow(mic),1],'.csv',sep=''))
    best = best[,c(1,2, ncol(best))]
    colnames(best) = c('V1','V2','Class')
    best$Class = factor(best$Class)

    worst= read.csv2(paste('../datasets/pca/microarray/',mic[1,1],'.csv',sep=''))
    worst = worst[,c(1,2, ncol(worst))]
    colnames(worst) = c('V1','V2','Class')
    worst$Class = factor(worst$Class)

    p1 = ggplot(best, aes(x=V1, y=V2, color=Class)) +
        geom_point() + 
        geom_rug() +
        scale_color_tableau(palette="colorblind10")+
        scale_fill_tableau(palette="colorblind10")+
        labs(title="Melhor dispers\u{00E3}o acumulada em duas componentes \n microarray",
             x=paste("Componente 1 (",round(porBest[1], 2),"% )"), y = paste("Componente 2 (",round(porBest[2],2),"% )"))+
        theme(plot.title = element_text(hjust = 0.5))

    p2 = ggplot(worst, aes(x=V1, y=V2, color=Class)) +
        geom_point() + geom_rug()+
        scale_color_tableau(palette="colorblind10")+
        scale_fill_tableau(palette="colorblind10")+
        labs(title="Pior dispers\u{00E3}o acumulada em duas componentes \n microarray",
             x=paste("Componente 1 (",round(porWorst[1], 2),"% )"), y = paste("Componente 2 (",round(porWorst[2],2),"% )")) + 
        theme(plot.title = element_text(hjust = 0.5))



    mi = cbind( name=as.character(df.pca.comp.mirnaseq[,1]), data.frame(sum=rowSums(df.pca.comp.mirnaseq[,c(2,3)])))
    mi = mi[order(mi[,2]),]
    porWorst = df.pca.comp.mirnaseq[df.pca.comp.mirnaseq[,1] == as.character(mi[1,1]),2:3] * 100
    porBest = df.pca.comp.mirnaseq[ df.pca.comp.mirnaseq[,1] == as.character(mi[nrow(mi),1]),2:3] * 100


    best= read.csv2(paste('../datasets/pca/miRNA-Seq/',mi[nrow(mi),1],'.csv',sep=''))
    best = best[,c(1,2, ncol(best))]
    colnames(best) = c('V1','V2','Class')
    best$Class = factor(best$Class)

    worst= read.csv2(paste('../datasets/pca/miRNA-Seq/',mi[1,1],'.csv',sep=''))
    worst = worst[,c(1,2, ncol(worst))]
    colnames(worst) = c('V1','V2','Class')
    worst$Class = factor(worst$Class)

    p3 = ggplot(best, aes(x=V1, y=V2, color=Class)) +
        geom_point() + 
        geom_rug()+
        scale_color_tableau(palette="colorblind10")+
        scale_fill_tableau(palette="colorblind10")+
        labs(title="Melhor dispers\u{00E3}o acumulada em duas componentes \n  miRNA-Seq",
             x=paste("Componente 1 (",round(porBest[1], 2),"% )"), y = paste("Componente 2 (",round(porBest[2],2),"% )"))+
        theme(plot.title = element_text(hjust = 0.5))

    p4 = ggplot(worst, aes(x=V1, y=V2, color=Class)) +
        geom_point() + geom_rug()+
        scale_color_tableau(palette="colorblind10")+
        scale_fill_tableau(palette="colorblind10")+
        labs(title="Pior dispers\u{00E3}o acumulada em duas componentes \n miRNA-Seq",
            x=paste("Componente 1 (",round(porWorst[1], 2),"% )"), y = paste("Componente 2 (",round(porWorst[2],2),"% )")) + 
        theme(plot.title = element_text(hjust = 0.5))




    rna = cbind( name=as.character(df.pca.comp.rnaseqv2[,1]), data.frame(sum=rowSums(df.pca.comp.rnaseqv2[,c(2,3)])))
    rna = rna[order(rna[,2]),]
    porWorst = df.pca.comp.rnaseqv2[df.pca.comp.rnaseqv2[,1] == as.character(rna[1,1]),2:3] * 100
    porBest = df.pca.comp.rnaseqv2[ df.pca.comp.rnaseqv2[,1] == as.character(rna[nrow(rna),1]),2:3] * 100


    best= read.csv2(paste('../datasets/pca/RNA-Seq/',rna[nrow(rna),1],'.csv',sep=''))
    best = best[,c(1,2, ncol(best))]
    colnames(best) = c('V1','V2','Class')
    best$Class = factor(best$Class)

    worst= read.csv2(paste('../datasets/pca/RNA-Seq/',rna[1,1],'.csv',sep=''))
    worst = worst[,c(1,2, ncol(worst))]
    colnames(worst) = c('V1','V2','Class')
    worst$Class = factor(worst$Class)

    p5 = ggplot(best, aes(x=V1, y=V2, color=Class)) +
        geom_point() + geom_rug()+
        scale_color_tableau(palette="colorblind10")+
        scale_fill_tableau(palette="colorblind10")+
        labs(title="Melhor dispers\u{00E3}o acumulada em duas componentes \n RNA-Seq",
            x=paste("Componente 1 (",round(porBest[1], 2),"% )"), y = paste("Componente 2 (",round(porBest[2],2),"% )"))+
        theme(plot.title = element_text(hjust = 0.5))



    p6 = ggplot(worst, aes(x=V1, y=V2, color=Class)) +
        geom_point() + geom_rug()+
        scale_color_tableau(palette="colorblind10")+
        scale_fill_tableau(palette="colorblind10")+
        labs(title="Pior dispers\u{00E3}o acumulada em duas componentes \n RNA-Seq",
            x=paste("Componente 1 (",round(porWorst[1], 2),"% )"), y = paste("Componente 2 (",round(porWorst[2],2),"% )")) + 
        theme(plot.title = element_text(hjust = 0.5))

    multiplot(p1,p3, p5, p2, p4, p6, cols=2, file="../output/analysis/pca.2comp.diff.pdf", save=T)
}

pca.all.diff = function(result){

    pcamean =  c(mean(result[result[,1]=='microarray','dim.PCA']),
                 mean(result[result[,1]=='mirnaseq','dim.PCA']),
                 mean(result[result[,1]=='rnaseqv2','dim.PCA']))
    allmean =  c(mean(result[result[,1]=='microarray','dim.Raw']),
                 mean(result[result[,1]=='mirnaseq','dim.Raw']),
                 mean(result[result[,1]=='rnaseqv2','dim.Raw']))

    ids = c('microarray', 'miRNA-Seq', 'RNA-Seq')
    pca.dt= data.frame(ids, pcamean)
    all.dt=data.frame(ids, allmean)

    p1 = ggplot(all.dt, aes(x=ids, y=allmean))+
        xlab("Conjunto de Dados") +
        ylab("N\u{00FA}mero m\u{00E9}dio de Atributos") +
        ggtitle("N\u{00FA}mero m\u{00E9}dio de atributos no conjunto de dados original ") +
        scale_color_tableau(palette="colorblind10")+
        scale_fill_tableau(palette="colorblind10")+
        geom_bar(stat = 'identity')+
        geom_text(stat='identity',
                  aes(label=format(allmean, digits=5, drop0trailing=TRUE), y= allmean ),vjust=-.5, fontface='bold')+
        theme(plot.title = element_text(hjust = 0.5, face='bold'), axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"), legend.text=element_text(size=14), legend.title=element_text(size=14, face='bold'))


    p2 = ggplot(pca.dt, aes(x=ids, y=pcamean))+
        xlab("Conjunto de Dados") +
        ylab("N\u{00FA}mero m\u{00E9}dio de Atributos (Componentes)") +
        scale_color_tableau(palette="colorblind10")+
        scale_fill_tableau(palette="colorblind10")+
        ggtitle("N\u{00FA}mero m\u{00E9}dio de atributos no conjunto de dados ap\u{00F3}s PCA-95% ") +
        geom_bar(stat = 'identity') +
        geom_text(stat='identity',
                  aes(label=format(pcamean, digits=4, drop0trailing=TRUE), y=pcamean ),vjust=-.5, fontface='bold') +
        theme(plot.title = element_text(hjust = 0.5, face='bold'), axis.text=element_text(size=14),
            axis.title=element_text(size=14,face="bold"), legend.text=element_text(size=14), legend.title=element_text(size=14, face='bold'))

    multiplot(p1, p2,cols=2)

    ggsave("../output/analysis/pca.all.diff.pdf")    

}

features.pca = function(result){
    levels(result$type) = c('RNA-Seq','miRNA-Seq','microarray')
    ggplot(result, aes(x=id, y=dim.PCA)) +
        xlab("Conjunto de Dados") +
        ylab("N\u{00FA}mero de Componentes") +
        ggtitle("Quantidade de componentes para 95% de dispers\u{00E3}o ") +
        geom_bar(stat = 'identity') + 
        scale_color_tableau(palette="colorblind10")+
        scale_fill_tableau(palette="colorblind10")+
        facet_grid(type ~ .)+
     theme(plot.title = element_text(hjust = 0.5, face='bold'), axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"), legend.text=element_text(size=14), legend.title=element_text(size=14, face='bold'))


    ggsave("../output/analysis/features.pca.pdf")    
}

disp.pca = function(){

    df.pca.comp.microarray = read.csv2('../output/percentage100PcaComponents.microarray.csv')
    df.pca.comp.mirnaseq = read.csv2('../output/percentage100PcaComponents.mirnaseq.csv')
    df.pca.comp.rnaseqv2 = read.csv2('../output/percentage100PcaComponents.rnaseqv2.csv')

    microarray = data.frame(dis=sort(apply(df.pca.comp.microarray[,c(2,3)],1, sum))*100)
    mirnaseq = data.frame(dis=sort(apply(df.pca.comp.mirnaseq[,c(2,3)],1,sum))*100)
    rnaseq = data.frame(dis=sort(apply(df.pca.comp.rnaseqv2[,c(2,3)],1,sum))*100)
    result= rbind(
        cbind(type=rep('microarray',dim(microarray)[1]), microarray),
        cbind(type=rep('miRNA-Seq', dim(mirnaseq)[1]),mirnaseq), 
        cbind(type=rep('RNA-Seq',dim(rnaseq)[1]),rnaseq))

   result[,'id'] = 1:86
   result$type = factor(result$type, levels=c('RNA-Seq','miRNA-Seq','microarray'))

    ggplot(result, aes(x=id, y=dis)) +
        xlab("Conjunto de Dados") +
        ylab("Contribui\u{00E7}\u{00E3}o na Dispers\u{00E3}o (%)") +
        ggtitle("Porcentagem de contribuicao na dispers\u{00E3}o \n nas duas primeiras componentes") +
        scale_color_tableau(palette="colorblind10")+
        scale_fill_tableau(palette="colorblind10")+
        geom_bar(stat = 'identity') + 
        facet_grid(type ~ .)+

        theme(plot.title = element_text(hjust = 0.5, face='bold'), axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"), legend.text=element_text(size=14), legend.title=element_text(size=14, face='bold'))



    ggsave("../output/analysis/disp.pca.pdf")    

}

agg.pca.all = function(){


    dataName = list.files(path = "../datasets/pca/", pattern = '*.csv', recursive=T)
    dataPath = paste('../datasets/pca/', dataName, sep='')
    dataPathraw = paste('../datasets/raw/', dataName, sep='')

    result = c()
    i=0
    for(f in dataPath){
        i=i+1
        data = read.csv2(f)
        dataraw = read.csv2(dataPathraw[i])
        result = rbind(result, data.frame(type=strsplit(strsplit(f,'/')[[1]][5],'\\.')[[1]][2], 
                                          dim.PCA=dim(data)[2]-1, dim.Raw=dim(dataraw)[2]-1))
        print(result[i,])
    }

    result = result[order(result$dim.PCA),]
    result = rbind(
                   result[result[,1]=='microarray',],
                   result[result[,1]=='mirnaseq',],
                   result[result[,1]=='rnaseqv2',])
    result[,'id'] = 1:86

    return(result)

    #grep('a', c('aa.c','c.d a','c','a'))
}

diff.classif = function(){
    
    metadata = read.table('../output/finalMetadata.dat')
    data = metadata[,1:11]
    data$data.name = c(rep('RNA-Seq',length(grep('rnaseqv2',data$data.name))), 
                       rep('miRNA-Seq',length(grep('mirnaseq',data$data.name))), 
                       rep('microarray',length(grep('microarray',data$data.name))))
    aux = data[,c('data.name','classif.svm.ber.test.mean', 'classif.randomForest.ber.test.mean', 'wilcox.test.p.value')]
    aux$diff = (1-aux$classif.svm.ber.test.mean) - (1-aux$classif.randomForest.ber.test.mean)
    aux = aux[order(aux[,c('diff')]),]
    aux$id = 1:86
    aux$test = ifelse(aux$wilcox.test.p.value < 0.05,'Dif. Estat\u00edstica', 'Sem Dif. Estat\u00edstica' )
    aux$test = factor(aux$test, levels=c('Sem Dif. Estat\u00edstica','Dif. Estat\u00edstica'))
    
     p1=ggplot(aux, aes(x=id, y=diff, fill=data.name,alpha=test)) +
        xlab("Conjunto de Dados") +
        ylab("Diferen\u{00E7}a BAC (SVM - RF)") +
        ggtitle("Diferen\u{00E7}a de desempenho (BAC) entre SVM e RF") +
        geom_bar(stat = 'identity')+
        scale_color_tableau(palette="colorblind10")+
        scale_fill_tableau(palette="colorblind10")+
        scale_alpha_discrete(range=c(0.4,1))+
        annotate("text", x = 75, y =0.35, label = "Ganho com SVM", size=5 , fontface="bold")+
        annotate("text", x = 15, y =-0.35, label = "Ganho com RF", size=5 , fontface="bold")+
        guides(fill=guide_legend(title="Tipo do Conjunto de Dados"), alpha=guide_legend(title='Transpar\u{00EA}ncia'))+
        theme(plot.title = element_text(hjust = 0.5, face='bold'), axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"), legend.text=element_text(size=14), legend.title=element_text(size=14, face='bold'))

    
    multiplot(p1, cols=1)

    ggsave("../output/analysis/diff.classif.pdf",plot=p1, width=14, height=7, dpi=300) 

}

table.datasets <- function(){


    dataName = list.files(path = "../datasets/raw/", pattern = '*.csv', recursive=T)
    dataPathraw = paste('../datasets/raw/', dataName, sep='')
   
    Tipo = c()
    Nome=c()
    ID=c()
    Qtd.Atributos=c()
    Qtd.Exemplos=c()
    Classes=c()

    for(f in dataPathraw){
        i=i+1
        dataraw = read.csv2(f)
        Tipo = c(Tipo, strsplit(f,'/')[[1]][4])
        Nome = c(Nome, strsplit(strsplit(f,'/')[[1]][5],'\\.')[[1]][1])
        ID = c(ID, paste('#',i,sep=' '))
        Qtd.Atributos = c(Qtd.Atributos, dim(dataraw)[1])
        Qtd.Exemplos = c( Qtd.Exemplos, dim(dataraw)[2])
        Classes = c( Classes, paste(sort(table(dataraw$class), decreasing=TRUE), "", collapse=""))
    }

    df = data.frame(ID=ID, Nome=Nome, Tipo=Tipo, Qtd.Atributos=Qtd.Atributos, Qtd.Exemplos=Qtd.Exemplos, Classes=Classes);
    xtable(df)
    return(df)
}

table.diff.SVMRF <- function(){

    casas=4
    metadata = read.table('../output/finalMetadata.dat')
    data = metadata[,1:11]
    data$data.name = c(rep('RNA-Seq',length(grep('rnaseqv2',data$data.name))), 
                       rep('miRNA-Seq',length(grep('mirnaseq',data$data.name))), 
                       rep('microarray',length(grep('microarray',data$data.name))))
    aux = data[,c('data.name','classif.svm.ber.test.mean', 'classif.randomForest.ber.test.mean', 'wilcox.test.p.value', 'classif.svm.timetrain.test.mean', 'classif.randomForest.timetrain.test.mean')]
    aux$SVM.BAC = round(1-aux$classif.svm.ber.test.mean,casas)
    aux$RF.BAC =round( 1-aux$classif.randomForest.ber.test.mean,casas)
    aux$Teste.Estatistico = ifelse(aux$wilcox.test.p.value > 0.05,'Sem Dif.', 'Dif.' )
    aux$ID = paste("#",1:86)
    aux$SVM.TEMPO = round(aux$classif.svm.timetrain.test.mean,casas)
    aux$RF.TEMPO = round(aux$classif.randomForest.timetrain.test.mean,casas)

    df=aux[,c('SVM.BAC','RF.BAC','SVM.TEMPO','RF.TEMPO','Teste.Estatistico')]
    colnames(df) = c("SVM BAC","RF BAC","SVM tempo (s)","RF tempo (s)","Teste Estatístico")
    rownames(df) = aux$ID
    xtable(df)
}

#    metadata = read.table('../output/finalMetadata.dat')

 #   types = c('^general.','^statistical.','^infotheo.','^model.based.', '^landmarking.', '^discriminant.')

  #  l = list()
  #  i=0
  #  for (type in types){
  #      i = i +1
  #      ids = grep(type,colnames(metadata))
  #      print(colnames(metadata[,ids]))
#
#        l[[i]] = colnames(metadata[,ids])
#    }



result = agg.pca.all()
#diff.classif()
#pca.all.diff(result)
#disp.pca()
#features.pca(result)
#pca.2comp.diff()

