This is a Matlab implementation of the method described in the following paper: R. J. Araújo, K. Fernandes, J. S. Cardoso, "Sparse Multibending Snakes", IEEE Transactions on Image Processing.

This implementation extends the code provided by Prof. Dirk-Jan Kroon in the Mathworks File Exchange, such that you need to download its contribution to the directory root:
https://www.mathworks.com/matlabcentral/fileexchange/28149-snake-active-contour

One of the examples provided (nr. 4) also requires a region growing implementation by Prof. Dirk-Jan
Kroon, which can be found at: 
https://www.mathworks.com/matlabcentral/fileexchange/19084-region-growing

The folders containing these pieces of code should be placed in the folder where this readme file is.

We implemented our optimization routine in C++, such that you have to use Matlab's mex function to compile the .cpp file into a binary MEX-file, callable from MATLAB.

This has been tested in Matlab R2015a and R2018a. For any problem or doubt do not hesitate in contacting me: rjtaraujo@gmail.com
