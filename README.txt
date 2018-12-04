This is a Matlab implementation of the method described in the following paper: R. J. Araújo, K. Fernandes, J. S. Cardoso, "Sparse Multibending Snakes", IEEE Transactions on Image Processing.

This implementation extends the code provided by Prof. Dirk-Jan Kroon in the Mathworks File Exchange, such that you need to download its contribution to the directory root:
https://www.mathworks.com/matlabcentral/fileexchange/28149-snake-active-contour

One of the examples provided also requires a region growing implementation by Prof. Dirk-Jan
Kroon, which can be found in: 
https://www.mathworks.com/matlabcentral/fileexchange/19084-region-growing

We provide our optimization routine implemented in C++, such that you have to use the mex
functionality of MATLAB in order to compile it. See more info at mex -setup

The folders containing these pieces of code should be placed in the folder where this readme file and the remaining files of our implementation are.

