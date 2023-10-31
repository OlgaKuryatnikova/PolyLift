These are implementations of the hierarchies from the paper, where all polynomials are coded via the function sdpvar in Yalmip. 
This approach is convenient for the Lasserre hierarchy and the sparse hierarchy from the paper. 
Schmudgen-based certificates can become slow in construction because of sdpvar polynomials (multiplying polynomials for the Schmudgen terms is the bottleneck). To speed up the construction (not the solution time, is stays more or less the same), we also code the polynomials using some functions from SOSTOOLS version 1.00 and DIGS instead of sdpvar. The resulting codes are in the folder Schmudgen_faster_codes. 

