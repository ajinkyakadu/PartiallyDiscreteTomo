Partially Discrete Tomography with Parametric Level-Set Method
==========================================

This toolbox solves the computational tomography problem where the object is partially discrete
in nature, that means - there is a variation part, called background, where properties varies between 0 and tau < 1, and an anomaly of constant value 1. This kind of objects occur in many applications, from image processing to medical imaging. We use parametric level-set method to represent the anomaly and the smooth regularization (i.e. tikhonov) for background. We compare the proposed method with Total-variation method, [DART](http://ieeexplore.ieee.org/document/5738333/)
and [P-DART](https://www.ncbi.nlm.nih.gov/pubmed/22417923).

Installation
------------

This toolbox requires Matlab version *R2016a or later*.  In
particular we require following packages to be pre-install to run this toolbox.

1. [ASTRA-Toolbox](https://github.com/astra-toolbox/astra-toolbox)  
2. [SPOT Operator](https://github.com/mpf/spot)  

Please run startup.m file before running the scripts in example folder.

Paper
------------

The theory and experiment details are available in the paper "A parametric level-set method for partially discrete tomography". The paper is available on ArXiv : [arXiv:1704.00568](https://arxiv.org/abs/1704.00568)


Feel free to contact me with any questions, suggestions, etc.

Ajinkya Kadu - [a.a.kadu@uu.nl](a.a.kadu@uu.nl)
