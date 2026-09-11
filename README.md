# PBSDC
Proximal bundle method for block-separable nonsmooth DC optimization

PBSDC [1] is a proximal double bundle solver (Fortran 95) for block-separable nonsmooth DC programming (difference of two convex functions) by Jenni Lampainen and Kaisa Joki. Within the code, we provide two variants of the proposed method. The first, PBSDC-crit, converges to an approximate $\varepsilon$-critical point, while the second, PBSDC-stat, converges to an approximate Clarke $\varepsilon$-stationary point.

During each round of 'main iteration' it is possible to utilize OpenMP to calculate subproblems in parallel. In order to do this, you need to use '-fopenmp' in Makefile.

The software utilizes code PVMM by Prof. Ladislav Luksan that is licensed by the GNU Lesser General Public License (LGPL). In addition, code PLQDF1 by Prof. Ladislav Luksan is used to solve quadratic direction finding problem.

The software also includes PBDC [2] and DBDC [3] methods by Kaisa Joki and both methods can be used in their original form in the code. Both methods are licensed by the MIT License.

The software is free for academic teaching and research purposes but I ask you to refer the reference given below if you use it. To use the software modify tpbsdc.f95 and functions.f95 as needed. If you have any questions conserning the software, please contact directly the author Jenni Lampainen (email: jmlamp@utu.fi).

## Code include:                                                                     
         
   tpbsdc.f95         - Main program for PBSDC 

   pbsdc.f95 	      - PBSDC method
   
   constants.f95      - Double precision (also some parameters) 
   
   bundle1.f95        - Bundle of DC component f_1
   
   bundle2.f95        - Bundle of DC component f_2                                    
        
   functions.f95      - User-specified DC components f_1 and f_2 together with subgradients of DC components. Contains also user-specified initial values for parameters.       
   
   fun.f95 	          - Defines objective funtion and gradient of the norm minimization problem

   norm_min.f95 	  - Solver for the norm minimization problem
   
   pbdc.f95 	      - PBDC method by K. Joki

   dbdc.f95 	      - DBDC method by K. Joki
	
   plqdf1.f 	      - Quadratic solver by L. Luksan
  
   pvmm.f             - Variable metric method by L. Luksan
   
   mqsubs.f 	      - Basic modules for PVMM by L. Luksan
   
   pqsubs.f 	      - Matrix modules for PVMM by L. Luksan
	
   Makefile 	      - Makefile

## Program usage: 
To use the code:
1. In the tpbsdc.f95 file, modify the parameters in lines 253–282. The most important parameters are solver_ind (1 = the new PBSDC method, 2 = the DBDC method, 3 = the PBDC method) and optimality_condition (1 = the PBSDC-crit variant, 2 = the PBSDC-stat variant). In addition, for example name of the starting point file and the output file can be specified later in tpbsdc.f95.
2. In functions.f95 file, define the DC functions in f1 and f2, and their subgradients in sugradient_f1 and subgradient_f2.
3. Run Makefile by typing "make".
4. Finally, just type "tpbsdc.exe".

The algorithm outputs a TXT file containing, for example, number of blocks, size of block, problem dimension, number of function evaluations, number of subgradient evaluations for the first and second DC components, CPU time, and function value.
   
## References:                                                                                         
Reference to PBSDC:

[1] J. Lampainen, K. Joki, A. M. Bagirov, S. Taheri and M. M. Mäkelä: "Proximal bundle method for block-separable nonsmooth DC optimization". Under review, (2026).

Reference to PBDC:

[2] Kaisa Joki, Adil M. Bagirov, Napsu Karmitsa and Marko M. Mäkelä: "A proximal bundle method for nonsmooth DC optimization utilizing nonconvex cutting planes". J. Glob. Optim. 68(3), 501-535, (2017).

Reference to DBDC:

[3] Kaisa Joki, Adil M. Bagirov, Napsu Karmitsa, Marko M. Mäkelä and Sona Taheri: "Double bundle method for finding Clarke stationary points in nonsmooth DC programming". SIAM J. Optim. 28(2), 1892-1919, (2018).
