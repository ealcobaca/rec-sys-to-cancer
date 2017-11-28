import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import os


def applyPCA(dirRead, dirSave):
    if not os.path.exists(dirSave):
        os.makedirs(dirSave)

    files = []
    filesS = []
    pcaVariance = []
    for file in os.listdir(dirRead):
        if file.endswith(".csv"):
            files.append(dirRead+'/'+ file)
            filesS.append(dirSave+'/' + file)
            pcaVariance.append([file.split(sep='.csv')[0]])
    i=-1
    for file, fileS in zip(files, filesS):
        i = i+1
        print('reading ... ['+str(file)+']')
        df = pd.read_csv(file, sep=';', decimal=",")

        y = df['Class']
        del df['Class']

        print('applying PCA ...')
        pca = PCA()
        pca.fit(df)
        counter = 0
        tot = 0
        for value in pca.explained_variance_ratio_:
            tot = tot + value
            counter = counter + 1
            if tot >= 0.95:
               break
        pcaVariance[i] = pcaVariance[i] + pca.explained_variance_ratio_[1:100].tolist()


        print('Number of components: ' + str(counter) + '')
        print('Variance percentage of components: '+ str(tot) + '')

        #print(df.shape)
        pca.n_components =  counter
        dfNew = pca.fit_transform(df)
        dfNew = pd.DataFrame(dfNew)
        #print(dfNew.shape)
        dfNew['Class'] = y

        print('saving ... ['+fileS+']\n\n')
        dfNew.to_csv(path_or_buf =fileS, sep=';', decimal=',', index=False)

    dfVariance = pd.DataFrame(pcaVariance)
    dfVariance.to_csv(
            path_or_buf='../output/percentage100PcaComponents.'+pcaVariance[0][0].split('.')[1]+'.csv',
            sep=';', decimal=',', index=False)


def main():
    dirsRead = ['../datasets/zscore/microarray', '../datasets/zscore/miRNA-Seq', '../datasets/zscore/RNA-Seq']
    dirsSave = ['../datasets/pca/microarray', '../datasets/pca/miRNA-Seq', '../datasets/pca/RNA-Seq']

    for i in range(0,len(dirsRead)):
        applyPCA(dirsRead[i],dirsSave[i])


if __name__ == "__main__":
    main()

