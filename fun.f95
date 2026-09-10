
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !| |                                                                                  | |       
        !| |   fun.f95 contains functions needed in the solver PVMM by Ladislav Luksan:       | |
        !| |                                                                                  | |    
        !| |   PVMM by Ladislav Luksan is needed in the norm minimization solver and it uses  | | 
        !| |   also PQSUBS and MQSUBS by Ladislav Luksan. PVMM is a VARIABLE METRIC ALGORITHM | |
        !| |   for UNCONSTRAINED and LINEARLY CONSTRAINED OPTIMIZATION.                       | |
        !| |                                                                                  | |
        !| |   The subroutine PVMM together with PQSUBS and MQSUBS are licensed by            | |
        !| |   the GNU Lesser General Public License (LGPL).                                  | |
        !| |                                                                                  | |
        !| |    --------------------------------------------------------------------------    | |
        !| |                                                                                  | |       
        !| |        FUNDER(NF, X, f_val, G) : Defines the objective function of the norm      | |
        !| |                                  minimization problem and its gradient           | |
        !| |                                                                                  | |
        !| |    --------------------------------------------------------------------------    | |
        !| |                                                                                  | |
        !| |                    **..**..* EMPTY SUBROUTINES *..**..**                         | |
        !| |                                                                                  | |
        !| |        FUN(NF,KA,X,FA)                                                           | |
        !| |        DFUN(NF,KA,X,GA)                                                          | |       
        !| |        OBJ(NF,X,FF)                                                              | |       
        !| |        DOBJ(NF,X,GF)                                                             | |       
        !| |        CON(NF,KC,X,FC)                                                           | |
        !| |        DCON(NF,KC,X,GC)                                                          | |       
        !| |                                                                                  | |
        !| |                                                                                  | |       
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*       

  
      SUBROUTINE FUNDER(NF, X, f_val, G)
      !
      ! Defines the objective function '|| f_matrix*X ||^2' of the norm minimization problem and its gradient '2 * g_matrix * X'. 
      ! 'f_matrix' contains subgradients from the bundles B_1 and B_2 (frst are the ones from B_1 and then from B_2).
      ! Needed in PVMM by Ladislav Luksan
      !
      ! INPUT:   * NF    : the number of variables in the norm minimization problem
      !          * X     : an estimate to the solution of the norm minimization problem
      !
      ! OUTPUT:  * f_val : the value of the objective function at the point 'X' in the norm minimization problem
      !          * G     : the subgradient of the objective function at the point 'X' in the norm minimization problem
      !
      ! NOTICE: The dimension of 'X' and 'G' is 'NF'.
      !
      !************************************************************************************
        USE constants, ONLY : dp                           ! double precision (i.e. accuracy)
        USE norm_min, ONLY  : n_orig, f_matrix, g_matrix   ! matrices which are used to define the objective function and its gradient 
                                                           ! 'n_orig' is the number of variables (i.e. dimension) in the original problem
      !************************************************************************************     
        IMPLICIT NONE
        
      ! ***** VARIABLES NEEDED FROM USER **************************************************
        INTEGER, INTENT(IN) :: NF                        ! the number of variables of the objective function in the norm minimization problem
        REAL(KIND=dp), INTENT(OUT) :: f_val              ! the value of the objective function at the point x
        REAL(KIND=dp), DIMENSION(NF), INTENT(IN)  :: X   ! an estimate to the solution
        REAL(KIND=dp), DIMENSION(NF), INTENT(OUT) :: G   ! the subgradient of the objective function at the point x
        
      ! ***** LOCAL VARIABLES **************************************************************    
        REAL(KIND=dp), DIMENSION(n_orig) :: vect_obj ! 'help' vector
        REAL(KIND=dp), DIMENSION(NF) :: vect_grad ! 'help' vector
        REAL(KIND=dp) :: help                     ! 'help' vector
        INTEGER :: i,j
        
        !-----------------------------------------------------------------------
        !*** Calulation of the objective function value f_val ***
        !-----------------------------------------------------------------------
        
         vect_obj = 0.0_dp  ! initialization
         
         ! Calculates: f_matrix * X
         DO j = 1, NF
            DO i = 1, n_orig
                 vect_obj(i) = vect_obj(i) + f_matrix(i,j)*X(j)
            END DO
         END DO 
         
         ! Objective: || f_matrix * X ||^2  
         f_val = DOT_PRODUCT(vect_obj, vect_obj)         

        !-----------------------------------------------------------------------
        !*** Calculation of the gradient G ***  
        !-----------------------------------------------------------------------
        
        vect_grad = 0.0_dp  ! initialization
        
        ! Calculates: g_matrix * X      
        DO j=1,NF
            DO i = 1, NF 
                vect_grad(i) = vect_grad(i) + g_matrix(i,j)*X(j)    
            END DO
        END DO
        
        ! Gradient: 2 * g_matrix * X 
        G = 2.0_dp * vect_grad                      


      END SUBROUTINE FUNDER
      
      
!************************************************************************
!       EMPTY SUBROUTINES (fortran 77 code)
!       Needed in PVMM by Ladislav Luksan

      SUBROUTINE FUN(NF,KA,X,FA)
      INTEGER KA,NF
      DOUBLE PRECISION FA,X(*)
      KA=NF
      FA=X(1)
      RETURN
      END
      
      SUBROUTINE DFUN(NF,KA,X,GA)
      INTEGER KA,NF
      DOUBLE PRECISION GA(*),X(*)
      KA=NF
      GA(1)=X(1)
      RETURN
      END
      
      SUBROUTINE OBJ(NF,X,FF)
      INTEGER NF
      DOUBLE PRECISION X(*),FF
      NF=1
      FF=X(1)
      RETURN
      END
      
      SUBROUTINE DOBJ(NF,X,GF)
      INTEGER NF
      DOUBLE PRECISION X(*),GF(*)
      NF=1
      GF(1)=X(1)
      RETURN
      END
      
      SUBROUTINE CON(NF,KC,X,FC)
      INTEGER NF,KC
      DOUBLE PRECISION X(*),FC
      KC=NF
      FC=X(1)
      RETURN
      END
      
      SUBROUTINE DCON(NF,KC,X,GC)
      INTEGER NF,KC
      DOUBLE PRECISION X(*),GC(*)
      KC=NF
      GC(1)=X(1)
      RETURN
      END
!************************************************************************     
  