all: prerpocessing base-level meta-level analysis

help:
	@echo ">> Please use \`make <target>' where <target> is one of <<"
	@echo "  preprocessing      to remove build files"
	@echo "  base-level         to run the program"
	@echo "  meta-level         to run the program"
	@echo "  analysis           to run the program"
	@echo "  all                to run the program"
	@echo "  claen              to run the program"
	@echo "  help               to run the program"

clean:
	@ echo ">> Cleaning files : '*~' '*.pyc' '*.pyo' <<"
	@ find . -name '*~' -exec rm --force  {} +
	@ find . -name '*.pyc' -exec rm --force {} +
	@ find . -name '*.pyo' -exec rm --force {} +

requirements:
	@ echo "Requirements for this application:\n"
	@ echo "Python 3.6"
	@ echo "    Pandas"
	@ echo "    Sklearn"
	@ echo "    Sklearn"
	@ echo "R 3.4.2"
	@ echo "    mlr"
	@ echo "    mfe"
	@ echo "    ggplot2"
	@ echo "    ggthemes"
	@ echo "    scales"

prerpocessing:
	@ echo ">> Running Preprocessing Steps  <<\n\n"
	@ echo ">> [Step 1] zscore transformation  <<\n\n"
	@ python zscore.py
	@ echo ">> [Step 2] principal component analysis PCA-95% <<\n\n"
	@ python fsPCA.py

base-level:
	@ echo ">> Running Base-Level Steps <<\n\n"
	@ echo ">> [Step 1] Extract Meta-Features  <<\n\n"
	@ R CMD BATCH getMfeFeatures.R
	@ echo ">> [Step 2] Running SVM and RF (base-level)  <<\n\n"
	@ R CMD BATCH runML.R
	@ echo ">> [Step 3] Building Meta-Base  <<\n\n"
	@ R CMD BATCH generetingMetaBase.R

meta-level:
	@ echo ">> Running Meta-Level Steps <<\n\n"
	@ echo ">> [Step 1] Building and Test Models  <<\n\n"
	@ R CMD BATCH runMLMetadata.R

analysis:
	@ echo ">> Data Analysis  <<\n\n"
	@ R CMD BATCH analysis.R

