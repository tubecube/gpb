#!--coding:utf-8--
import sklearn.metrics as metric
import matplotlib.pyplot as plt

path='/Users/tubecube/Documents/graduate'

def plot_roc(y_test, y_prob, title, is_save=False):
	"""
	绘制ROC曲线图
	:param 默认不保存
	:return: 
	"""
	fpr, tpr, thresholds = metric.roc_curve(y_test, y_prob, pos_label=1)
	auc = metric.roc_auc_score(y_test, y_prob)

	plt.figure()
	lw = 2
	plt.plot(fpr, tpr, color='darkorange', lw=lw, label='ROC curve (area = %0.2f)' % auc)
	plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
	plt.xlim([0.0, 1.0])
	plt.ylim([0.0, 1.05])
	plt.xlabel('False Positive Rate')
	plt.ylabel('True Positive Rate')

	plt.legend(loc="best")
	plt.title('AUC_{0}: auc={1:0.2f}'.format(title, auc))
	plt.show()
	if is_save:
		plt.savefig(path+'/AUC_{0}_{1:0.2f}.jpg'.format(title, auc))

def plot_prc(y_test, y_prob, is_save=False):
	"""
	适用于样本不平衡分类文同的模型评估
	:param  默认不保存
	:return: 
	"""
	if self.y_score:
		average_precision = metric.average_precision_score(self.y_test, self.y_score)
		precision, recall, _ = metric.precision_recall_curve(self.y_test, self.y_score)
	else:
		average_precision = metric.average_precision_score(self.y_test, self.y_prob[:, 1])
		precision, recall, _ = metric.precision_recall_curve(self.y_test, self.y_prob[:,1])

	plt.step(recall, precision, color='b', alpha=0.2, where='post')
	plt.fill_between(recall, precision, step='post', alpha=0.2, color='b')

	plt.xlabel('Recall')
	plt.ylabel('Precision')
	plt.ylim([0.0, 1.05])
	plt.xlim([0.0, 1.0])
	plt.title('PRC_{0}_{1}: AP={2:0.2f}'.format(self.model_name, self.param, average_precision))
	plt.show()

	if is_save:
		plt.savefig(self.path+'/PRC_{0}_{1}_{2:0.2f}.jpg'.format(self.model_name, self.param, average_precision))

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
	filename = sys.argv[1]
	labels, probs = read_label_prob(filename)
	plot_roc(labels, probs, filename)
