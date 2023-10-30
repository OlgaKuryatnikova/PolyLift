These are implementations of the hierarchies from the paper, where all polynomials are coded via the function sdpvar in Yalmip. 
This approach is fast and convenient for the Lasserre hierarchy and the sparse hierarchy from the paper. 
Schmudgen-based certificates are slow-to-construct when using sdpvar polynomials (obtaining the Schmudgen terms is the bottleneck). To speed up the construction (not the solution time, is stays more or less the same), we code the polynomials using some functions from SOSTOOLS version 1.00 (first reference below) and DIGS (second reference below) instead of sdpvar. The resulting codes are in the folder Schmudgen_faster_codes. 

References:
1. A. Papachristodoulou, J. Anderson, G. Valmorbida, S. Prajna, P. Seiler, P. A. Parrilo, M. M. Peet and D. Jagt, SOSTOOLS: Sum of squares optimization toolbox for MATLAB, http://arxiv.org/abs/1310.4716, 2021, available from https://github.com/oxfordcontrol/SOSTOOLS.
2. B. Ghaddar, J.C. Vera, & M.F. Anjos, A dynamic inequality generation scheme for polynomial programming. Math. Program. 156, 21–57 (2016). https://doi.org/10.1007/s10107-015-0870-9.

