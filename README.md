# Joint Demosaicing and Denoising of Noisy Bayer Images with ADMM

``Hanlin Tan, Xiangrong Zeng, Shiming Lai,  Yu Liu and Maojun Zhang``

``College of Information System and Management, National University of Defense Technology, China``


This project is a demo for our ICIP 2017 paper [Joint Demosaicing and Denoising of noisy Bayer Images with ADMM](https://www.researchgate.net/publication/317058420_JOINT_DEMOSAICING_AND_DENOISING_OF_NOISY_BAYER_IMAGES_WITH_ADMM?_iepl%5BviewId%5D=nNdQs0DA16qpA3QVVefsanKG&_iepl%5BprofilePublicationItemVariant%5D=default&_iepl%5Bcontexts%5D%5B0%5D=prfpi&_iepl%5BtargetEntityId%5D=PB%3A317058420&_iepl%5BinteractionType%5D=publicationTitle). The project compare three algorithms, DeepJoint[\[2\]][2], FlexISP[\[3\]][3] and the proposed ADMM algorithm, on two datasets: Kodak and McMaster[\[1\]][1].

Note that we slightly modified I/O code (mainly redirecting input and output paths) of FlexISP and DeepJoint(Demosaicnet) to perform comparisons. However, the kernel part of those algorithms, including all parameters, remains the same as their authors provided.

## Folder Structure

| Folder Name | Explanation |
|:-:|:-:|
| data | Datasets Kodak and McMater[\[1\]][1] are here. |
| demo\_joint\_isp | Code of the proposed ADMM method. |
| demosaicnet | DeepJoint[\[2\]][2] code. |
| flexisp\_demosaicking\_demo | FlexISP[\[3\]][3] code. |
| res | Results. |
| utils | Code needed for PSNR computation and Bayer mask computation. |
| wrappers | Wrappers for comparative algorithms |



## Usage

In summary, useful entrance files are:

| File Name | Explanation |
|:-:|:-:|
| `compareAlgorithms.m` | The entrance file. Run this file with Matlab and the results are stored in `compareResults.mat` in `res/{dataset}` folder. |
| `analyze.m` | Produce figures and tables of results. This file is called automatically by `compareAlgorithms.m`. Results will be saved in `res/{dataset}` folder. |
| `combineResultImages.m` | Reproduce figures of visual results using data in `res` folder. |



Simply run `compareAlgorithms.m`, which will generate a result mat file in `res/{dataset}` folder.

+ You can control the behavior of the code by modifying structure variable `conf`.

+ Note `DeepJoint` method depends on `caffe`. You can turn off evaluation of that algorithm by setting `conf.enableDeepJoint` to `false` in the control section of `compareAlgorithms.m` if you do not have caffe installed.

```matlab
%% control
conf.debug = false;
conf.enableADMM = true;      % enable admm method
conf.enableFlexISP = true;   % enable flex isp method
conf.enableDeepJoint = false; % enable deep joint method
conf.cpuorgpu = 'gpu';    % 'gpu' or 'cpu', wheter to use gpu in deep joint
conf.tmpdir = '~/tmp/';   % a directory to store temporary output
```

It will call `analyze(dataset, sigma, nAlgos)` to produce result images, where `dataset` and `sigma` are in accordance with those in `compareAlgorithms.m`. For exmaple,

```
analyze('Kodak', 0, 3)
```

Optionally, run `combineResultImages(dataset)` to combine result images of different noise levels into single view, which we used in our paper. For exmaple,

```
combineResultImages('Kodak')
```


## References
If you use the code in your work, please cite our paper:

```
@inproceedings{tan2017joint, 
  title={Joint Demosaicing and Denoising of noisy Bayer Images with ADMM},
  author={Tan, Hanlin and Zeng, Xiangrong and Lai, Shiming and Liu, Yu and Zhang, Maojun},
  booktitle={Image Processing (ICIP), 2017 IEEE International Conference on},
  pages={2951--2955},
  year={2017},
  organization={IEEE}
}
```

The comparison code directly use the dataset and code from the following papers:

[\[1\] Zhang, Lei, et al. "Color demosaicking by local directional interpolation and nonlocal adaptive thresholding." Journal of Electronic imaging 20.2 (2011): 023016-023016.][1]

[\[2\] Gharbi, Michaël, et al. "Deep joint demosaicking and denoising." ACM Transactions on Graphics (TOG) 35.6 (2016): 191.][2]

[\[3\] Heide, Felix, et al. "Flexisp: A flexible camera image processing framework." ACM Transactions on Graphics (TOG) 33.6 (2014): 231.][3]

[1]:http://ira.lib.polyu.edu.hk/bitstream/10397/6039/1/Zhang_Color_Demosaicking_Local.pdf
[2]:http://dl.acm.org/citation.cfm?id=2982399
[3]:http://dl.acm.org/citation.cfm?id=2661260