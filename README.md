# PolyLift
Matlab code for the project https://arxiv.org/abs/1909.06689. 

**Examples:** examples used in the paper
**Yalmip_codes**: codes for all certificates from the paper;  Schmudgen-based certificates are somewhat slower in construction (not in solution) than Schmudgen_faster_codes.
**Schmudgen_faster_codes:** codes for Schmudgen-based certificates with faster construction times (due to faster polynomial multplication). Schmudgen_faster_codes use some functions from SOSTOOLS version 1.00 and DIGS.
**Support_functions:** support functions needed to construct and run the certificates. 

Besides the above codes, the certificates use Yalmip, Mosek (for SDP), Gurobi (for LP and SOCP).  

References:
1.  SOSTOOLS: A. Papachristodoulou, J. Anderson, G. Valmorbida, S. Prajna, P. Seiler, P. A. Parrilo, M. M. Peet and D. Jagt, SOSTOOLS: Sum of squares optimization toolbox for MATLAB, http://arxiv.org/abs/1310.4716, 2021, available from https://github.com/oxfordcontrol/SOSTOOLS.
2. DIGS: B. Ghaddar, J.C. Vera, & M.F. Anjos, A dynamic inequality generation scheme for polynomial programming. Math. Program. 156, 21–57 (2016). https://doi.org/10.1007/s10107-015-0870-9.
3.  Yalmip: https://yalmip.github.io/reference/lofberg2004/
4.  Mosek: https://www.mosek.com/
5.  Gurobi: https://www.gurobi.com/
