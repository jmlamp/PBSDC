      MODULE constants
        IMPLICIT NONE

        ! ** Double precision (i.e accuracy) **
        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12)  
        
        ! ** Integer variables that can be used if in the solver PVMM of the norm minimization problem we want to know some extra information ** 
		!     - 'NRES'  number of RESTARTS.
		!     - 'NDEC'  number of MATRIX DECOMPOSITION.
		!     - 'NREM'  number of CONSTRAINT DELETIONS.
		!     - 'NADD'  number of CONSTRAINT ADDITIONS.
		!     - 'NIT'  number of ITERATIONS.
		!     - 'NFV'  number of FUNCTION EVALUATIONS.
		!     - 'NFG'  number of GRADIENT EVALUATIONS.
		!     - 'NFH'  number of HESSIAN EVALUATIONS.		
		
        !INTEGER, SAVE :: NRES,NDEC,NREM,NADD,NIT,NFV,NFG,NFH       
  
      END MODULE constants