# Stereo
这里的代码都是与stereo相关的代码，包括自己制作的数据集，以及练习所写的代码。

## matlab-test-dataset

这个文件夹下含有的是我利用matlab所合成的两个平面，一个是在前的**斜面**，另一个是背景的**平面**，含有**相机的内外参**，以及**groundtruth**数据。

编写的语言自然是matlab，使用的版本是Mac上的matlab 2015。

##matlab-l2norm-minimun surface

这个代码对应的是**Efficient Minimal-Surface Regularizationof Perspective Depth Maps in Variational Stereo**，这篇文章采用的是最小化曲面面积实现的透视模型的正则项，此处我使用的是仅仅是简单的l2-norm,而非作者使用的l2，1-norm，主要是因为l2-norm 实现简单，可以使用梯度下降法去求解，而不是采用复杂的primal-dual，可以帮助去测试代码中别的部分是否正确，之所以将其改写为matlab，是因为matlab调试的时候实在是太轻松了。
l0-norm是想使用零范数去实现的失败版本，l2-norm6等也是同样。

##matlab-minimum surface

这个代码对应的是**Efficient Minimal-Surface Regularizationof Perspective Depth Maps in Variational Stereo**仅仅是实现了用matlab改写的过程，后续会使用自己的数据集进行测试。

## python-l0smooth
这个代码是对应**Image Smoothing via L0 Gradient Minimization**，其中原作者给出的代码是matlab所写的，这里因为需要就将其改写成了python版本。

其中使用的opencv库都是直接用anaconda管理的，下载更新都十分的方便。

