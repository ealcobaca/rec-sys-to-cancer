import pandas as pd
from sklearn import preprocessing
import glob, os


def normDatasets(dirRead, dirSave):
    if not os.path.exists(dirSave):
        os.makedirs(dirSave)

    files = []
    filesS = []
    for file in os.listdir(dirRead):
        if file.endswith(".csv"):
            files.append(dirRead+'/'+ file)
            filesS.append(dirSave+'/' + file)

    for file, fileS in zip(files, filesS):
        # print(file)
        # print(fileS)
        print('reading ... [' + file + ']')
        df = pd.read_csv(file, sep=';', decimal=",")

        y = df['class']
        #y = np.split(y, y.shape[0])
        del df['class']

        X = df.values

        #preprocessecing
        standardized_X = preprocessing.scale(X)

        dfNew = pd.DataFrame(standardized_X)
        dfNew['Class'] = y

        print('saving ... ['+fileS+']\n\n')
        dfNew.to_csv(fileS, sep=';', decimal=',', index=False)


def main():
    dirsRead = ['../datasets/raw/microarray', '../datasets/raw/miRNA-Seq', '../datasets/raw/RNA-Seq']
    dirsSave = ['../datasets/zscore/microarray', '../datasets/zscore/miRNA-Seq', '../datasets/zscore/RNA-Seq']

    for i in range(0,len(dirsRead)):
        normDatasets(dirsRead[i],dirsSave[i])

if __name__ == "__main__":
    main()

