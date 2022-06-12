import os
os.chdir('D:/E/博士/R_程序/HGSOC/DataGEOtogather')


##### 2022.2.16 ######
# https://machinelearningmastery.com/rfe-feature-selection-in-python/
# check scikit-learn version
import sklearn
print(sklearn.__version__)


## RFE 是一种变换。要使用它，首先为类配置通过“ estimator ”参数指定的所选算法和通过“ n_features_to_select ”参数选择的特征数量。


import numpy as np
from pandas import DataFrame
from sklearn.feature_selection import RFE

import pandas as pd
# load data
dataframe = pd.read_table("matrix_DEtrain.txt")
# sample name
samplenames = np.array(dataframe.columns)
# featire name
featurenames = np.array(dataframe.index)
# data form
esetarray = np.array(dataframe)
esetarray = esetarray.transpose()

# regression
#sampletype = [0]*10+[1]*10
#sampletype0 =  pd.read_table("trainlab.txt")
#sampletype = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
#                       1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
#                       1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
#                       1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
#                       1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
#                       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
#                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
#                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
#                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
#                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
#                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
#                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

sampletype = [0]*8+[1]*6+[0]*14+[1]*5+[0]*1

p = 2022

#############################   SVM-RFE    ###############################

from sklearn import svm
clf = svm.SVC(kernel='linear', random_state=p)  

# =============================================================================
# rfe = RFE(estimator = clf,           # 基分类器
#           n_features_to_select = 1, # 选择特征个数
#           step = 1,                  # 每次迭代移除的特征个数 
#           verbose = 1                # 显示中间过程
#           )
# =============================================================================

rfe = RFE(clf, n_features_to_select = 1)
rfe.fit(esetarray, sampletype)
rfe.score(esetarray, sampletype)
# 和传参对应，所选择的属性的个数
print(rfe.n_features_)
# 打印的是相应位置上属性的排名
print(rfe.ranking_)
# 属性选择的一种模糊表示，选择的是true，未选择的是false
print(rfe.support_)
# 第1个属相的排名
print(rfe.ranking_[1])
# 外部估计函数的相关信息
print(rfe.estimator_)


print("Features sorted by their rank:")
print(sorted(zip(map(lambda x: round(x, 4), rfe.ranking_), featurenames)))
result = sorted(zip(map(lambda x: round(x, 4), rfe.ranking_), featurenames))
resultframe = DataFrame(result)
resultframe.to_csv("DataPython/ranklist_SVMrfe.txt", sep="\t")


#############################  DecisionTree-RFE    ###############################

from sklearn.tree import DecisionTreeClassifier
# define the method
rfe = RFE(estimator=DecisionTreeClassifier(random_state=p), n_features_to_select=1)
# fit the model
rfe.fit(esetarray, sampletype)
# 类拟合后，可以通过“ support_ ”属性查看输入变量的选择，该属性为每个输入变量提供True或False。
rfe.support_
# 和传参对应，所选择的属性的个数
print(rfe.n_features_)
# 打印的是相应位置上属性的排名
print(rfe.ranking_)
# 属性选择的一种模糊表示，选择的是true，未选择的是false
print(rfe.support_)
# 第1个属相的排名
print(rfe.ranking_[1])
# 外部估计函数的相关信息
print(rfe.estimator_)


print("Features sorted by their rank:")
print(sorted(zip(map(lambda x: round(x, 4), rfe.ranking_), featurenames)))
result = sorted(zip(map(lambda x: round(x, 4), rfe.ranking_), featurenames))
resultframe = DataFrame(result)
resultframe.to_csv("DataPython/ranklist_DTrfe.txt", sep="\t")


#############################  RandomForestClassifier-RFE    ###############################

from sklearn.ensemble import RandomForestClassifier
rfe = RFE(estimator=RandomForestClassifier(random_state=p), n_features_to_select=1)
# fit the model
rfe.fit(esetarray, sampletype)
rfe.score(esetarray, sampletype)
# 类拟合后，可以通过“ support_ ”属性查看输入变量的选择，该属性为每个输入变量提供True或False。
rfe.support_
# 和传参对应，所选择的属性的个数
print(rfe.n_features_)
# 打印的是相应位置上属性的排名
print(rfe.ranking_)
# 属性选择的一种模糊表示，选择的是true，未选择的是false
print(rfe.support_)
# 第1个属相的排名
print(rfe.ranking_[1])
# 外部估计函数的相关信息
print(rfe.estimator_)


print("Features sorted by their rank:")
print(sorted(zip(map(lambda x: round(x, 4), rfe.ranking_), featurenames)))
result = sorted(zip(map(lambda x: round(x, 4), rfe.ranking_), featurenames))
resultframe = DataFrame(result)
resultframe.to_csv("DataPython/ranklist_RFrfe.txt", sep="\t")


#############################  GradientBoostingClassifier-RFE    ###############################

from sklearn.ensemble import GradientBoostingClassifier
rfe = RFE(estimator=GradientBoostingClassifier(random_state=p), n_features_to_select=1)
# fit the model
rfe.fit(esetarray, sampletype)
rfe.score(esetarray, sampletype)
# 类拟合后，可以通过“ support_ ”属性查看输入变量的选择，该属性为每个输入变量提供True或False。
rfe.support_
# 和传参对应，所选择的属性的个数
print(rfe.n_features_)
# 打印的是相应位置上属性的排名
print(rfe.ranking_)
# 属性选择的一种模糊表示，选择的是true，未选择的是false
print(rfe.support_)
# 第1个属相的排名
print(rfe.ranking_[1])
# 外部估计函数的相关信息
print(rfe.estimator_)



# now print out the features in order of ranking
from operator import itemgetter
features = featurenames
for x, y in (sorted(zip(rfe.ranking_ , features), key=itemgetter(0))):
    print(x, y)


# =============================================================================
# # ok, this time let's choose the top 10 featues and use them for the model
# n_features_to_select = 10
# rfe = RFE(GradientBoostingClassifier, n_features_to_select=n_features_to_select)
# rfe.fit(esetarray, sampletype)
# # use the model to predict the prices for the test data
# predictions = rfe.predict(X_test)
# # write out CSV submission file
# output = pd.DataFrame({"Id":test_data.index, target:predictions})
# output.to_csv('submission.csv', index=False)
# 
# =============================================================================


print("Features sorted by their rank:")
print(sorted(zip(map(lambda x: round(x, 4), rfe.ranking_), featurenames)))
result = sorted(zip(map(lambda x: round(x, 4), rfe.ranking_), featurenames))
resultframe = DataFrame(result)
resultframe.to_csv("DataPython/ranklist_GBMrfe.txt", sep="\t")


#############################  AdaBoost-RFE    ###############################
from sklearn.ensemble import AdaBoostRegressor as AdaBoost
clf = AdaBoost(random_state=p)  
rfe = RFE(clf, n_features_to_select=1) 
# fit the model
rfe.fit(esetarray, sampletype)
rfe.score(esetarray, sampletype)
# 类拟合后，可以通过“ support_ ”属性查看输入变量的选择，该属性为每个输入变量提供True或False。
rfe.support_
# 和传参对应，所选择的属性的个数
print(rfe.n_features_)
# 打印的是相应位置上属性的排名
print(rfe.ranking_)
# 属性选择的一种模糊表示，选择的是true，未选择的是false
print(rfe.support_)
# 第1个属相的排名
print(rfe.ranking_[1])
# 外部估计函数的相关信息
print(rfe.estimator_)



# now print out the features in order of ranking
from operator import itemgetter
features = featurenames
for x, y in (sorted(zip(rfe.ranking_ , features), key=itemgetter(0))):
    print(x, y)

print("Features sorted by their rank:")
print(sorted(zip(map(lambda x: round(x, 4), rfe.ranking_), featurenames)))
result = sorted(zip(map(lambda x: round(x, 4), rfe.ranking_), featurenames))
resultframe = DataFrame(result)
resultframe.to_csv("DataPython/ranklist_ABrfe.txt", sep="\t")


#############################  xgboost-RFE    ###############################

import xgboost as xgb

# =============================================================================
# https://stackoverflow.com/questions/55658234/python-3-modulenotfounderror-no-module-named-xgboost
# import sys
# !{sys.executable} -m pip install xgboost
# 
# =============================================================================

clf = xgb.XGBRegressor(random_state=p)   
rfe = RFE(clf, n_features_to_select=1) 
# fit the model
rfe.fit(esetarray, sampletype)
rfe.score(esetarray, sampletype)
# 类拟合后，可以通过“ support_ ”属性查看输入变量的选择，该属性为每个输入变量提供True或False。
rfe.support_
# 和传参对应，所选择的属性的个数
print(rfe.n_features_)
# 打印的是相应位置上属性的排名
print(rfe.ranking_)
# 属性选择的一种模糊表示，选择的是true，未选择的是false
print(rfe.support_)
# 第1个属相的排名
print(rfe.ranking_[1])
# 外部估计函数的相关信息
print(rfe.estimator_)



# now print out the features in order of ranking
from operator import itemgetter
features = featurenames
for x, y in (sorted(zip(rfe.ranking_ , features), key=itemgetter(0))):
    print(x, y)

print("Features sorted by their rank:")
print(sorted(zip(map(lambda x: round(x, 4), rfe.ranking_), featurenames)))
result = sorted(zip(map(lambda x: round(x, 4), rfe.ranking_), featurenames))
resultframe = DataFrame(result)
resultframe.to_csv("DataPython/ranklist_XGBrfe.txt", sep="\t")


#############################  NaiveBayes-RFE    ###############################
#from sklearn.naive_bayes import GaussianNB
#clf = GaussianNB()    
from sklearn.naive_bayes import BernoulliNB
clf = BernoulliNB()   
rfe = RFE(clf, n_features_to_select=1) 
# fit the model
rfe.fit(esetarray, sampletype)
rfe.score(esetarray, sampletype)
# 类拟合后，可以通过“ support_ ”属性查看输入变量的选择，该属性为每个输入变量提供True或False。
rfe.support_
# 和传参对应，所选择的属性的个数
print(rfe.n_features_)
# 打印的是相应位置上属性的排名
print(rfe.ranking_)
# 属性选择的一种模糊表示，选择的是true，未选择的是false
print(rfe.support_)
# 第1个属相的排名
print(rfe.ranking_[1])
# 外部估计函数的相关信息
print(rfe.estimator_)


print("Features sorted by their rank:")
print(sorted(zip(map(lambda x: round(x, 4), rfe.ranking_), featurenames)))
result = sorted(zip(map(lambda x: round(x, 4), rfe.ranking_), featurenames))
resultframe = DataFrame(result)
resultframe.to_csv("DataPython/ranklist_NBrfe.txt", sep="\t")






#############################  KNN-RFE    ###############################
#from sklearn.neighbors import KNeighborsClassifier
#clf = KNeighborsClassifier(n_neighbors = 5)   
#rfe = RFE(clf, n_features_to_select=1) 
## fit the model
#rfe.fit(esetarray, sampletype)
#rfe.score(esetarray, sampletype)
## 类拟合后，可以通过“ support_ ”属性查看输入变量的选择，该属性为每个输入变量提供True或False。
#rfe.support_
## 和传参对应，所选择的属性的个数
#print(rfe.n_features_)
## 打印的是相应位置上属性的排名
#print(rfe.ranking_)
## 属性选择的一种模糊表示，选择的是true，未选择的是false
#print(rfe.support_)
## 第1个属相的排名
#print(rfe.ranking_[1])
## 外部估计函数的相关信息
#print(rfe.estimator_)
#
#
#
## now print out the features in order of ranking
#from operator import itemgetter
#features = featurenames
#for x, y in (sorted(zip(rfe.ranking_ , features), key=itemgetter(0))):
#    print(x, y)
#
#print("Features sorted by their rank:")
#print(sorted(zip(map(lambda x: round(x, 4), rfe.ranking_), featurenames)))
#result = sorted(zip(map(lambda x: round(x, 4), rfe.ranking_), featurenames))
#resultframe = DataFrame(result)
##resultframe.to_csv("D:/E/博士/R_程序/HGSOC/DataPython/ranklist_KNNbrfe.txt", sep="\t")










## explore the algorithm wrapped by RFE
#from numpy import mean
#from numpy import std
#from sklearn.datasets import make_classification
#from sklearn.model_selection import cross_val_score
#from sklearn.model_selection import RepeatedStratifiedKFold
#from sklearn.feature_selection import RFE
#from sklearn.linear_model import LogisticRegression
#from sklearn.linear_model import Perceptron
#from sklearn.tree import DecisionTreeClassifier
#from sklearn.ensemble import RandomForestClassifier
#from sklearn.ensemble import GradientBoostingClassifier
#from sklearn.pipeline import Pipeline
#from matplotlib import pyplot
#
## get the dataset
#def get_dataset():
#	X, y = make_classification(n_samples=1000, n_features=10, n_informative=5, n_redundant=5, random_state=1)
#	return X, y
#
## get a list of models to evaluate
#def get_models():
#	models = dict()
#	# lr
#	rfe = RFE(estimator=LogisticRegression(), n_features_to_select=5)
#	model = DecisionTreeClassifier()
#	models['lr'] = Pipeline(steps=[('s',rfe),('m',model)])
#	# perceptron
#	rfe = RFE(estimator=Perceptron(), n_features_to_select=5)
#	model = DecisionTreeClassifier()
#	models['per'] = Pipeline(steps=[('s',rfe),('m',model)])
#	# cart
#	rfe = RFE(estimator=DecisionTreeClassifier(), n_features_to_select=5)
#	model = DecisionTreeClassifier()
#	models['cart'] = Pipeline(steps=[('s',rfe),('m',model)])
#	# rf
#	rfe = RFE(estimator=RandomForestClassifier(), n_features_to_select=5)
#	model = DecisionTreeClassifier()
#	models['rf'] = Pipeline(steps=[('s',rfe),('m',model)])
#	# gbm
#	rfe = RFE(estimator=GradientBoostingClassifier(), n_features_to_select=5)
#	model = DecisionTreeClassifier()
#	models['gbm'] = Pipeline(steps=[('s',rfe),('m',model)])
#	return models
#
## evaluate a give model using cross-validation
#def evaluate_model(model, X, y):
#	cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=1, random_state=1)
#	scores = cross_val_score(model, X, y, scoring='accuracy', cv=cv, n_jobs=None)    # n_jobs=-1   并行 'n_jobs=None'
#	return scores
#
## define dataset
#X, y = get_dataset()
## get the models to evaluate
#models = get_models()
## evaluate the models and store results
#results, names = list(), list()
#for name, model in models.items():
#	scores = evaluate_model(model,esetarray, sampletype)
#	results.append(scores)
#	names.append(name)
#	print('>%s %.3f (%.3f)' % (name, mean(scores), std(scores)))
## plot model performance for comparison
#pyplot.boxplot(results, labels=names, showmeans=True)
#pyplot.show()






