#!--coding:utf-8--
from scipy import interp
import numpy as np
import sklearn.metrics as metric
import matplotlib.pyplot as plt

path='/Users/tubecube/Documents/graduate'

mean_tpr = np.zeros(100)
mean_fpr = np.linspace(0, 1, 100)
mean_auc = 0.0

def plot_mean(fpr, tpr, auc, title="Protein230(K=7)"):
	plt.figure()
	lw = 1
	plt.plot(fpr, tpr, color='darkorange', lw=lw, label='ROC curve (area = %0.2f)' % auc)
	plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
	plt.xlim([-0.05, 1.0])
	plt.ylim([0.0, 1.05])
	plt.xlabel('False Positive Rate')
	plt.ylabel('True Positive Rate')

	plt.legend(loc="best")
	plt.title('{0}: auc={1:0.4f}'.format(title, auc))
	plt.show()


def plot_roc(y_test, y_prob):
	"""
	绘制ROC曲线图
	:param 默认不保存
	:return: 
	"""
	global mean_tpr, mean_auc
	fpr, tpr, thresholds = metric.roc_curve(y_test, y_prob, pos_label=1,
	drop_intermediate=True)
	auc = metric.roc_auc_score(y_test, y_prob)
	mean_tpr += interp(mean_fpr, fpr, tpr)
	mean_auc += auc

def plot_prc(y_test, y_prob, title):
	"""
	适用于样本不平衡分类文同的模型评估
	:param  默认不保存
	:return: 
	"""
	average_precision = metric.average_precision_score(y_test, y_prob)
	precision, recall, _ = metric.precision_recall_curve(y_test, y_prob)

	plt.subplot(1,2,1)
	plt.step(recall, precision, color='b', alpha=0.2, where='post')
	# plt.fill_between(recall, precision, step='post', alpha=0.2, color='b')

	plt.xlabel('Recall')
	plt.ylabel('Precision')
	plt.xlim([-0.05, 1.0])
	plt.ylim([0.0, 1.05])
	plt.title('PRC_{0}: AP={1:0.2f}'.format(title, average_precision))
	plt.subplot(1,2,2)
	plt.plot(recall, precision)
	plt.show()

def read_label_prob(filename):
	labels = []
	probs = []
	for line in open(filename):
		a, b, label, prob = line.strip().split('\t')
		labels.append(int(label))
		probs.append(float(prob))
	return labels, probs

if __name__ == "__main__":
	import sys
	for filename in sys.argv[1:]:
		labels, probs = read_label_prob(filename)
		plot_roc(labels, probs)
	mean_tpr /= len(sys.argv[1:])
	mean_auc /= len(sys.argv[1:])
	mean_tpr[0] = 0
	plot_mean(mean_fpr, mean_tpr, mean_auc)
	# plot_prc(labels, probs, filename)
