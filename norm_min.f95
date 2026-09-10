      MODULE norm_min

        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*!
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |!
        !| |                                                                                  | |!
        !| |                  THE SOLVER FOR THE NORM MINIMIZATION PROBLEM                    | |! 
        !| |                                                                                  | |!
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |!
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*!   
        !|                                                                                      |!
        !|                                                                                      |!
        !|    Utilizes PVMM by Ladislav Luksan as a norm minimization solver. This subroutine   |! 
        !|    uses PQSUBS and MQSUBS by Ladislav Luksan. PVMM is a VARIABLE METRIC ALGORITHM    |!
        !|    for UNCONSTRAINED and LINEARLY CONSTRAINED OPTIMIZATION.                          |!
        !|                                                                                      |!
        !|    The subroutine PVMM together with PQSUBS and MQSUBS are licensed by               |!
        !|    the GNU Lesser General Public License (LGPL).                                     |!
        !|                                                                                      |!
        !|                                                                                      |!
        !|                                                                                      |!
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |!
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*!
        
        
        USE constants   ! Double precision (i.e. accuracy) and if needed the varibles: NRES,NDEC,NREM,NADD,NIT,NFV,NFG,NFH
        USE bundle1     ! The BUNDLE of the DC component f_1
        USE bundle2     ! The BUNDLE of the DC component f_2
                
        IMPLICIT NONE 
        
        EXTERNAL PVMM       ! Contains the Variable metric method PVVML by Ladislav Luksan and it is needed to solve 
                            ! the norm minimization problem at Step 3 of the 'main iteration' algorithm.
        
        EXTERNAL PQSUBS     ! Basic modules for PVVM (by Ladislav Luksan)
        EXTERNAL MQSUBS     ! Matrix modules for PVVM (by Ladislav Luksan) 
        
        INTEGER, SAVE :: n_orig                                      ! The number of variables (i.e. dimension) in the original problem
                                                                     ! and it is needed/used in the EXTERNAL subroutine FUNDER.

        ! Matrices which are used to define the objective function of the norm minimization problem and its gradient. 
        ! Both needed/used in the EXTERNAL subroutine FUNDER.
        REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE, SAVE :: f_matrix    ! used to define the objective function in the norm minimization problem.
        REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE, SAVE :: g_matrix   ! used to define gradient of the objective function in the norm minimization problem.
    
        
        CONTAINS
        
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..** |
        !| |                                                                                             | |
        !| |                            CONTAINS SUBROUTINES:                                            | | 
        !| |                                                                                             | |
        !| |    SOLVER FOR THE NORM PROBLEM             : norm_solver(bxi, obj, NF, n1, n2, B1, B2)      | | 
        !| |                                                                                             | |
        !| |    FORMS MATRICES USED IN THE NORM PROBLEM : norm_matrix(dimensio, n1, n2, B1, B2)          | |
        !| |                                                                                             | |       
        !| |                                                                                             | |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..** |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..



        !***************************************************************************************
        !  ----------------------------------------------------------------------------------  |
        !  |                                                                                |  |
        !  |                   SOLVER FOR THE NORM MINIMIZATION PROBLEM                     |  |
        !  |                                                                                |  |
        !  ----------------------------------------------------------------------------------  |
        !***************************************************************************************        
        
            SUBROUTINE norm_solver(bxi, obj, NF, n1, n2, B1, B2, user_n)
                !
                ! Solves the norm minimization problem at Step 3 of the 'main iteration' algorithm.
                !
                ! INPUT:  * 'B1' and 'B2':  The bundles B_1 and B_2 of the DC component f_1 and f_2 
                !         * 'n1' and 'n2':  The size of the bundles B_1 and B_2 (the 'current bundle element' is also taken into account in these values but 'aggregate element' is NOT!)
                !         * 'NF'         :  The number of variables in the objective function of the norm minimization problem (NOTICE: NF = n1 + n2)
                !         * 'user_n'     :  The dimension of the original problem
                !
                ! OUTPUT: * 'bxi':  The difference of the combined subgradients bxi = bxi_1^* - bxi_2^* = f_matrix * lambda ('lambda' is the approximate minimum to the norm minimization problem)
                !         * 'obj':  The value of the objective funtion at the solution 'lambda' (i.e. '|| bxi ||')
                !
            !***********************************************************************************
                IMPLICIT NONE
                
            !**************************** NEEDED FROM USER *************************************    
                TYPE(kimppu1), INTENT(IN) :: B1  ! the bundle B_1 for the DC component f_1
                TYPE(kimppu2), INTENT(IN) :: B2  ! the bundle B_2 for the DC component f_2
                                
                REAL(KIND=dp), DIMENSION(give_n_b1(B1)), INTENT(OUT) :: bxi     ! Output: bxi = bxi_1^* - bxi_2^* = f_matrix * lambda
                REAL(KIND=dp), INTENT(OUT) :: obj                               ! Output: the value of the objective function at the solution 'lambda' 
                
                INTEGER, INTENT(IN) :: NF        ! the number of variables in objective function of the norm min problem (NOTICE: NF = n1 + n2)
                INTEGER, INTENT(IN) :: user_n    ! the dimension of the original problem
                INTEGER, INTENT(IN) :: n1, n2    ! the size of the bundle B_1 and B_2 (the 'current bundle element' is also taken into account in these values but 'aggregate element' is NOT!)
                
            
            !***************************** LOCAL VARIABLES WHEN CALLING PVMML ******************

                ! .. Help variables ..                        
                INTEGER :: i,j,ind
                
                ! .. Other variables ..
                
                INTEGER, PARAMETER :: NB = 3  ! the simple bounds are accepted because NB > 0 
                INTEGER, PARAMETER :: NC = 2  ! the number of linear constraints (2 equality constraints)
                
                INTEGER :: IPRNT               ! Specifies print in PVMML:  IPRNT = 0  print is  suppressed
                                               !                            IPRNT = 1  basic print of final results
                                               !                            IPRNT = -1 extended print of final results 
                                               !                            IPRNT = 2  basic print of intermediate and final results
                                               !                            IPRNT = -2  extended print of intermediate and final results
                                               
                INTEGER :: ITERM               ! OUTPUT from PVMML: indicates the cause of termination              
                
                REAL(KIND=dp) :: GMAX, F                 ! OUTPUT from PVMML : F    - the value of the objective at the solution
                                                         !                     GMAX - maximum absolute value of a partial derivative of the objective function              
                
                INTEGER, DIMENSION(7) :: IPAR            ! INTEGER parameters (If IPAR = 0, then default values are used in PVMML)                
                
                REAL(KIND=dp), DIMENSION(NF) :: lambda   ! when INPUT:  a vector with the initial estimate to the solution
                                                         ! when OUTPUT: the approximation to the minimum               
                                                         
                REAL(KIND=dp), DIMENSION(7) :: RPAR      ! REAL parameters (IF RPAR = 0.0_dp, then default values are used in PVMML)


                ! .. Array Arguments in PVMML..
                INTEGER, DIMENSION(NF)       :: IX     ! the vector containing types of simple bounds (IX=1 because XL(I) <= X(I) for each variable)
                REAL(KIND=dp), DIMENSION(NF) :: XL     ! the vector containing lower bounds (XL(I) = 0.0_dp for each variable)
                REAL(KIND=dp), DIMENSION(NF) :: XU     ! vector containing upper bounds     (This is not needed)
                
                REAL(KIND=dp), DIMENSION(NC) :: CF     ! optional/auxiliary: the vector containing values of constraint functions
                INTEGER, DIMENSION(NC)       :: IC     ! the vector containing constraint types: 5 - the equality constraint (NOTICE: CL(I) = C(I) = CU(I) )
                REAL(KIND=dp), DIMENSION(NC) :: CL     ! the vector containing lower bounds for constraints 
                REAL(KIND=dp), DIMENSION(NC) :: CU     ! the vector containing upper bounds for constraints
                REAL(KIND=dp), DIMENSION(NF*NC) :: CG  ! the vector whose columns are normals of linear constraints                 

                
                ! .. Initialization ..

                n_orig = user_n    ! The dimension of the original problem is stored to the extra variable 'n_orig'
                
                ALLOCATE(f_matrix(user_n,NF),g_matrix(NF,NF))   ! Allocation of the correct sizes for 'f_matrix' and 'g_matrix'


                IPAR = 0         ! default values are used
                RPAR = 0.0_dp    ! default values are used
                IPRNT = 0        ! print (because IPRNT=0, print is suppressed)
                
                IX = 1            ! the simple bounds are lower bounds
                XL = 0.0_dp       ! values of the lower bounds
                XU = 10.0_dp      ! values of the upper bounds (NOT needed here!)
                
                !CF = (/  /)                 ! optional         
                IC = (/ 5, 5 /)              ! constraints are equality constraints
                CL = (/ 1.0_dp, 1.0_dp /)    ! lower bounds for constraints
                CU = (/ 1.0_dp, 1.0_dp /)    ! upper bounds for constraints
                
                CG = 0.0_dp         ! Initialization of CG          
                DO i = 1, n1        ! the first constraint (indices: 1,...,n1+n2)
                    CG(i) = 1.0_dp              
                END DO
                
                ind = 2*n1 + n2 + 1  
                DO i = ind, (ind + n2-1)  ! the second constraint (indices: n1+n2+1,...,2*n1+2*n2)
                    CG(i) = 1.0_dp
                END DO
                
                ! initial estimate to the solution
                lambda = 0.0_dp         ! initialization of lambda
                lambda(1) = 1.0_dp
                lambda(n1+1) = 1.0_dp
                
                !i=give_solution_ind(B2)
                !lambda(n1+i) = 1.0_dp
                            
                CALL norm_matrix(user_n, n1, n2, B1, B2) ! Calculates the matrices used to define the objective function and 
                                                         ! its gradient in the norm minimization problem
                
                ! Execution of the solver PVMML by Ladislav Luksan
                CALL PVMML(NF,NB,NC,lambda,IX,XL,XU,CF,IC,CL,CU,CG,IPAR,&     ! objective F = || f_matrix * lambda ||^2
                            & RPAR, F, GMAX, IPRNT, ITERM )
          

                IF ( F < 0.0_dp) THEN
                    obj = 0.0_dp
                ELSE
                    obj = SQRT(F)    ! the original objective function 'obj' = || f_matrix * lambda || 
                END IF

                bxi = 0.0_dp         ! initialization 
                
                DO j = 1, NF
                    DO i = 1, user_n
                       bxi(i) = bxi(i) + f_matrix(i,j)*lambda(j)
                    END DO
                END DO

                DEALLOCATE(f_matrix,g_matrix)   ! Deallocation of 'f_matrix' and 'g_matrix'           
                
         
            END SUBROUTINE norm_solver
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -       
            
            
            
        !***************************************************************************************
        !  ----------------------------------------------------------------------------------  |
        !  |                                                                                |  |
        !  |            CALCULATE MATRICES f_matrix and g_matrix USED IN FUNDER             |  |
        !  |                                                                                |  |
        !  ----------------------------------------------------------------------------------  |
        !***************************************************************************************
        
            SUBROUTINE norm_matrix(user_n, n1, n2, B1, B2) 
             ! Calculates the matrices used to define the objective of norm minimization problem
             ! and its gradient. These matrices are saved to the variables 'f_matrix' and 'g_matrix'
             !
             ! NOTICE: * 'user_n'       : the length of the subgradient
             !         * 'n1' and 'n2'  : the size of the bundles B_1 and B_2 (the 'current bundle element' is also taken into account in these values but 'aggregate element' is NOT!)
             !         * 'B1' and 'B2*  : the bundles B_1 and B_2 of the DC components f_1 and f_2
             !
                 INTEGER, INTENT(IN) :: user_n, n1, n2    ! user_n  - the length of the subgradients (i.e. dimension of the original problem)
                                                          ! n1 and n2 - the size of the bundles B_1 and B_2 (the 'current bundle element' is also taken into account in these values but 'aggregate element' is NOT!)
                   
                 INTEGER :: i, j, ind
     
                 ! .. Help variables ..
                 REAL(KIND=dp), DIMENSION(user_n) :: grad_b1, grad_b2                  
                 REAL(KIND=dp), DIMENSION(n1+n2,n1+n2) :: M             ! M = grad^T * grad  
                 REAL(KIND=dp), DIMENSION(user_n, n1+n2) :: grad        ! the matrix whose columns are subgradients
                                                                        ! contains all the subgradients (both from B_1 and B_2, but 'aggregate elements' are NOT included)
                                                                        ! first there are the ones from B_1 and then the ones from B_2
                 ! .. Bundles ..                                                        
                 TYPE(kimppu1), INTENT(IN) :: B1          ! the bundle of the DC component f_1
                 TYPE(kimppu2), INTENT(IN) :: B2          ! the bundle of the DC component f_2                                                                      

 
                 ! .. Matrix 'grad' is formed ..
                 DO i = 0, n1-1                         ! all subgradients from the bundle B_1 are looked through (except the 'aggregate element')
                    grad_b1 = give_subgrad_b1(B1,i)     ! the subgradient of the bundle element i from the bundle B_1
                    DO j = 1, user_n                    ! inserts the current grad_b1 into the matrix 'grad' (the starting position is (i+1):th column)
                       grad(j, i+1) = grad_b1(j)  
                    END DO                    
                 END DO 
               
                 ind = n1 + 1                           ! the starting position for subgradients from the bundle B_2
               
                 DO i = 0, n2-1                         ! all subgradients from the bundle B_2 are looked through (except the 'aggregate element')
                    grad_b2 = give_subgrad_b2(B2,i)     ! the subgradient of the bundle element i from the bundle B_2
                    DO j = 1, user_n                    ! inserts the current grad_b2 into the matrix 'grad' (NOTICE: the negative sign)
                       grad(j, ind + i) = - grad_b2(j)  ! (the position is (ind+i):th column in the matrix)
                    END DO                
                 END DO         
                 
                 ! .. Matrix 'grad' is saved to 'f_matrix' ..
                 DO j = 1, n1+n2
                    DO i = 1, user_n
                      f_matrix(i,j) = grad(i,j)
                    END DO
                 END DO              
                 
                 M =  MATMUL( TRANSPOSE(grad), grad )      ! M = grad^T * grad
                 
                 ! .. Matrix 'M' is saved to 'g_matrix' ..
                 DO j = 1, n1+n2
                    DO i = 1, n1+n2
                      g_matrix(i,j) = M(i,j)
                    END DO
                 END DO
       
             
            END SUBROUTINE norm_matrix      
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            
        

      END MODULE norm_min        