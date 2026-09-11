        MODULE functions    
      
        USE constants, ONLY   : dp   ! double precision (i.e. accuracy)    
        IMPLICIT NONE
        
      
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !| |                                                                                  | |
        !| |                                                                                  | |
        !| |                 INFORMATION SUPPLIED BY THE USER:                                | | 
        !| |                                                                                  | |
        !| |    * PROBLEM specification GIVEN by user in tpbsdc.f95:                          | |
        !| |                                                                                  | |       
        !| |        - The DC component f1:                    'problem1'                      | |
        !| |        - The DC component f2:                    'problem2'                      | |
        !| |        - the number of variables:                'user_n'                        | |  
        !| |        - the number of blocks                    'number_of_blocks'              | |       
        !| |        - the number of variables in one block    'size_of_block'                 | |
        !| |                                                                                  | |
        !| |    * Different PARAMETERS:                                                       | |
        !| |                                                                                  | |       
        !| |        - the stopping tolerance:      'user_crit_tol'                            | |
        !| |                                                                                  | |       
        !| |        GLOBAL PARAMETERS:                                                        | |       
        !| |        - the number of blocks                     'number_of_blocks'             | |       
        !| |        - the number of variables in one block     'size_of_block'                | |       
        !| |                                                                                  | |       
        !| |        MAIN ITERATION:                                                           | |       
        !| |        - the size of bundle B_1:      'user_size_b1'                             | |       
        !| |        - the size of bundle B_2:      'user_size_b2'                             | |       
        !| |        - the descent parameter:       'user_m'                                   | |       
        !| |        - the decrease parameter:      'user_c'                                   | |       
        !| |        - the decrease parameter:      'user_r_dec'                               | |       
        !| |        - the increase parameter:      'user_r_inc'                               | |           
        !| |        - the enlargement parameter:   'user_eps_1'                               | |       
        !| |                                                                                  | |       
        !| |        CLARKE STATIONARY ALGORITHM:                                              | |       
        !| |        - the size of bundle B:        'user_size'                                | | 
        !| |        - the proximity measure:       'user_eps'                                 | |               
        !| |        - the descent parameter:       'user_m_clarke'                            | |       
        !| |                                                                                  | |       
        !| |                                                                                  | |       
        !| |    * Computation of the value of the DC functions f_1 and f_2:                   | |
        !| |        - f1(y, problem1, user_n)   the value of DC component f_1 at a point y    | |
        !| |        - f2(y, problem2, user_n)   the value of DC component f_2 at a point y    | |           
        !| |                                                                                  | |               
        !| |                                                                                  | |               
        !| |    * Computation of the subgradient of the DC components f_1 and f_2:            | |
        !| |        - subgradient_f1(y, problem1, user_n)  the subgradient of f_1 at y        | |
        !| |        - subgradient_f2(y, problem2, user_n)  the subgradient of f_2 at y        | |       
        !| |                                                                                  | |
        !| |                                                                                  | |
        !| |                                                                                  | |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*

        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*|
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**.. |
        !| |                                                                                           | |
        !| |                                                                                           | |
        !| |                              SUBROUTINES AND FUNCTIONS:                                   | |       
        !| |                                                                                           | |       
        !| |    * Problem data allocation and deallocation:                                            | |
		!| |       - allocate_prob_data(nblock, sizeblock)                                             | |       
        !| |       - deallocate_prob_data()                                                            | |             
        !| |                                                                                           | |   
        !| |    * Calculates the sum used in f_1 and f_2 for parameter t_j:                            | |
        !| |       - summa(y, j, user_n)                                                               | |   
        !| |                                                                                           | |   
        !| |    * Computation of the value of the DC functions f_1 and f_2:                            | |
        !| |        - f1(y, problem1, user_n)   the value of DC component f_1 at a point y             | |
        !| |        - f2(y, problem2, user_n)   the value of DC component f_2 at a point y             | |	
        !| |                                                                                           | |		
        !| |    * Computation of the value of the block-separable DC functions f_1 and f_2:            | |
        !| |        - f1_block(y, problem1, user_n_block)   the value of DC component f_1 at a point y | |
        !| |        - f2_block(y, problem2, user_n_block)   the value of DC component f_2 at a point y | |                         
        !| |                                                                                           | |  
		!| |    * Computation of the subgradient of the DC components f_1 and f_2:                     | |
        !| |        - subgradient_f1(y, problem1, user_n)  the subgradient of f_1 at y                 | |
        !| |        - subgradient_f2(y, problem2, user_n)  the subgradient of f_2 at y                 | |
        !| |                                                                                           | |		
        !| |    * Computation of the subgradient of the block-separable DC components f_1 and f_2:     | |
        !| |        - subgradient_f1_block(y, problem1, user_n_block)  the subgradient of f_1 at y     | |
        !| |        - subgradient_f2_block(y, problem2, user_n_block)  the subgradient of f_2 at y     | |       
        !| |                                                                                           | |
        !| |                                                                                           | |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**.. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*|

		!*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !| |                                                                              | |
        !| |     To EXECUTE the bundle algorithm PBSDC the USER needs to DETERMINE:       | | 
        !| |                                                                              | |       
        !| |        f1(y, problem1, user_n)     - value of DC component f_1 at a point y  | |
        !| |        f2(y, problem2, user_n)     - value of DC component f_2 at a point y  | |
        !| |                                                                              | |
        !| |        subgradient_f1(y, problem1, user_n)   - subgradient of f_1 at y       | |
        !| |        subgradient_f2(y, problem2, user_n)   - subgradient of f_2 at y       | |
        !| |                                                                              | |
        !| |                                                                              | |       
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*

        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !| |                                                                                                              | |
        !| |                                    TEST PROBLEMS FROM PBSDC PAPER [1]:                                       | |
        !| |                                                                                                              | |
        !| |                                                                                                              | |
        !| |    Separable problems:                                                                                       | |
        !| |                                                                                                              | |
        !| |                 DC components f1 and f2                                                                      | |
        !| |               ---------------------------                                                                    | |       
        !| |                                                                                                              | |
        !| |    Problem 1:             S1, S2                                                                             | |
        !| |    Problem 2:             S1, S3                                                                             | |
        !| |    Problem 3:             S1, S4                                                                             | |
        !| |    Problem 4:             S2, S4                                                                             | |
        !| |    Problem 5:             S5, S3                                                                             | |
        !| |    Problem 6:             S5, S4                                                                             | |
        !| |    Problem 7:             S6, S4                                                                             | |
        !| |    Problem 8:             S7, S1                                                                             | |
        !| |    Problem 9:             S7, S4                                                                             | |
        !| |    Problem 10:            S7, S8                                                                             | |
        !| |                                                                                                              | |
        !| |                                                                                                              | |
        !| |    Block-separable problems:                                                                                 | |
        !| |                                                                                                              | |
        !| |                 DC components f1 and f2                                                                      | |
        !| |               ---------------------------                                                                    | |   
        !| |                                                                                                              | |
        !| |    Problem 11:            B1, B3                                                                             | |
        !| |    Problem 12:            B1, B5                                                                             | |
        !| |    Problem 13:            B1, B6                                                                             | |
        !| |    Problem 14:            B1, B7                                                                             | |
        !| |    Problem 15:            B1, B12                                                                            | |
        !| |    Problem 16:            B2, B5                                                                             | |
        !| |    Problem 17:            B2, B7                                                                             | |
        !| |    Problem 18:            B3, B2                                                                             | |
        !| |    Problem 19:            B3, B4                                                                             | |
        !| |    Problem 20:            B3, B5                                                                             | |
        !| |    Problem 21:            B3, B6                                                                             | |
        !| |    Problem 22:            B4, B5                                                                             | |
        !| |    Problem 23:            B4, B7                                                                             | |
        !| |    Problem 24:            B5, B7                                                                             | |
        !| |    Problem 25:            B8, B6                                                                             | |
        !| |    Problem 26:            B8, B7                                                                             | |
        !| |    Problem 27:            B8, B13                                                                            | |
        !| |    Problem 28:            B9, B7                                                                             | |
        !| |    Problem 29:            B9, B12                                                                            | |
        !| |    Problem 30:            B9, B13                                                                            | |
        !| |    Problem 31:            B10, B7                                                                            | |
        !| |    Problem 32:            B10, B12                                                                           | |
        !| |    Problem 33:            B10, B13                                                                           | |
        !| |    Problem 34:            B11, B3                                                                            | |
        !| |    Problem 35:            B11, B6                                                                            | |
        !| |    Problem 36:            B11, B8                                                                            | |
        !| |    Problem 37:            B11, B12                                                                           | |
        !| |    Problem 38:            B11, B13                                                                           | |
        !| |    Problem 39:            B11, B15                                                                           | |
        !| |    Problem 40:            B11, B16                                                                           | |
        !| |    Problem 41:            B14, B1                                                                            | |
        !| |    Problem 42:            B14, B2                                                                            | |
        !| |    Problem 43:            B14, B3                                                                            | |
        !| |    Problem 44:            B14, B4                                                                            | |
        !| |    Problem 45:            B14, B5                                                                            | |
        !| |    Problem 46:            B14, B6                                                                            | |
        !| |    Problem 47:            B14, B7                                                                            | |
        !| |    Problem 48:            B14, B8                                                                            | |
        !| |    Problem 49:            B14, B9                                                                            | |
        !| |    Problem 50:            B14, B12                                                                           | |
        !| |    Problem 51:            B14, B13                                                                           | |
        !| |    Problem 52:            B14, B16                                                                           | |
        !| |    Problem 53:            B15, B2                                                                            | |
        !| |    Problem 54:            B15, B4                                                                            | |
        !| |    Problem 55:            B15, B6                                                                            | |
        !| |    Problem 56:            B15, B13                                                                           | |
        !| |    Problem 57:            B15, B16                                                                           | |
        !| |    Problem 58:            B16, B6                                                                            | |
        !| |                                                                                                              | |
        !| |                                                                                                              | |
        !| |  [1] J. Lampainen, K. Joki, A. M. Bagirov, S. Taheri and M. M. Mäkelä: "Proximal bundle method for           | |
        !| |      block-separable nonsmooth DC optimization". Under review, (2026).                                       | |
        !| |                                                                                                              | |                                                                                                              | |
        !| |                                                                                                              | |                                                                                                              | |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        
        
        !--------------------------------------------------------------------------------
        ! ----------------------------------------------------------------------------- |
        ! |                  INFORMATION ABOUT PARAMETERS:                            | |
        ! ----------------------------------------------------------------------------- |
        !--------------------------------------------------------------------------------               

        !****************** GLOBAL PARAMETERS *******************************************
      
        INTEGER, SAVE :: number_of_blocks        ! The number of blocks in block-separable functions
        
        INTEGER, ALLOCATABLE, SAVE :: size_of_block(:)           ! The number of variables in one block

        !-------------------------------------------------------------------------------------------------      
        !__________________________________________________________________________________________
        !>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>>
        !****************** PARAMETRES NEEDED ONLY IN MAIN ITERATION AlGORITHM ********************
        
        INTEGER, PARAMETER :: user_size_b1 = -3                   ! The biggest possible size of the bundle B_1
                                                                  ! If user_size_b1 <= 0 then DEFAULT value MIN(user_n+5,100) is used 
                                                                   
        INTEGER, PARAMETER :: user_size_b2 = 3                    ! The biggest possible size of the bundle B_2
                                                                  ! If user_size_b2 <= 0 then DEFAULT value 3 is used  
              
        REAL(KIND=dp), PARAMETER :: user_m = 0.2_dp               ! The descent parameter r in PBSDC
                                                                  ! If user_m <= 0.0_dp .OR. user_m >= 1.0_dp then DEFAULT value 0.2_dp is used

        REAL(KIND=dp), PARAMETER :: user_c = -0.1_dp              ! The decrease parameter c2 in PBSDC
                                                                  ! If user_c <= 0.0_dp or user_c > 1.0_dp then DEFAULT value 0.1_dp is used
                                                                 
        REAL(KIND=dp), PARAMETER :: user_r_dec = -0.99_dp         ! The decrease parameter c1 in PBSDC
        
        !If user_r_dec <= 0.0_dp .OR. user_r_dec >= 1.0_dp then DEFAULT value is used.
        !                               
        !   DEFAULT value:                          
        !     If user_n < 10:           user_r_dec = 0.75_dp    
        !     If 10 <= user_n < 300:    user_r_dec = the first two decimals of n/(n+5)
        !     If user_n >= 300:         user_r_dec = 0.99_dp
        !
        !   Some examples of the DEFAULT value of the parameter 'user_r_dec':
        !     If user_n=10:     user_r_dec = 0.66_dp                          
        !     If user_n=20:     user_r_dec = 0.80_dp                         
        !     If user_n=25:     user_r_dec = 0.83_dp                         
        !     If user_n=50:     user_r_dec = 0.90_dp                         
        !     If user_n=100:    user_r_dec = 0.95_dp                         
        !     If user_n=150:    user_r_dec = 0.96_dp     
        !     If user_n=200:    user_r_dec = 0.97_dp                      
        !     If user_n=250:    user_r_dec = 0.98_dp    
        !
        
        REAL(KIND=dp), PARAMETER :: user_r_inc = (10.0_dp)**(7)       ! The increase parameter C in PBSDC
                                                                      ! If user_r_inc <= 1.0_dp then DEFAULT value (10.0_dp)**7 is used
 
        REAL(KIND=dp), PARAMETER :: user_eps_1 = 5*(10.0_dp)**(-5)    ! The enlargement parameter in PBSDC
                                                                      ! If user_eps_1 <= 0.0_dp .OR. user_eps_1 > 1.0_dp then DEFAULT value 5*(10.0_dp)**(-5) is used

 
        !____________________________________________________________________________________________                                                       
        !>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<
        !****************** PARAMETRES NEEDED ONLY IN CLARKE STATIONARY AlGORITHM *******************
        
        
         INTEGER, PARAMETER :: user_size = - 2               ! The biggest possible bundle size for the bundle used in 'Clarke stationary' algorithm
                                                             ! In this version of code user_n cannot be used here. If user_size < 2 then DEFAULT MIN(user_n*2,200) is used

         REAL(KIND=dp), PARAMETER :: user_m_clarke = 0.01_dp        ! The descent parameter: If user_m_clarke <= 0.0_dp .OR. user_m_clarke >= 1.0_dp 
                                                                    !                        then DEFAULT value 0.01_dp is used
                                                                    
          
         REAL(KIND=dp), PARAMETER :: user_eps = -(10.0_dp)**(-7)    ! The proximity measure: If user_eps <= 0.0_dp 
                                                                    !              then for PBSDCstat DEFAULT value (10.0_dp)**(-6) is used when n <= 50
                                                                    !                                               (10.0_dp)**(-5) is used when n > 50
                                                                    !              then for PBSDCcrit DEFAULT value 0.1_dp is used
          
         ! The (overestimated) Lipschitz constant of the DC component f_1 on the set {x\in\R^n | d(x, F_0) <= user_eps }      
        REAL(KIND=dp), PARAMETER :: user_L1 = 1000.0_dp
        
        ! The (overestimated) Lipschitz constant of the DC component f_2 on the set {x\in\R^n | d(x, F_0) <= user_eps }     
        REAL(KIND=dp), PARAMETER :: user_L2 = 1000.0_dp
        ! NOTICE: If user_L1 <= 0.0_dp or user_L2 <= 0.0_dp then DEFAULT value 1000.0_dp is used    

          
        !________________________________________________________________________________
        !>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<
        !****************** PARAMETER NEEDED IN STOPPING CONDITIONS ********************
              
        ! ** The stopping tolerance **      
        REAL(KIND=dp), PARAMETER :: user_crit_tol = -(10.0_dp)**(-5)   ! For PBSDCstat
                                                                       !  If user_crit_tol <= 0.0_dp then DEFAULT value (10.0_dp)**(-5) is used when n <=200
                                                                       !                                                (10.0_dp)**(-4) is used when n > 200
                                                                       ! For PBSDCcrit
                                                                       !  If n <= 0.0_dp then DEFAULT value is used
                                                                       !   DEFAULT value:  
                                                                       !     If user_n < 150:           user_crit_tol = user_n*0.005_dp      
                                                                       !     If 150 <= user_n <= 200:   user_crit_tol = user_n*0.015_dp      
                                                                       !     If user_n > 200:           user_crit_tol = user_n*0.05_dp 
        
        
        
        !________________________________________________________________________________
        !>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<
        
        CONTAINS

          !!----------------------------------------------------------------------------------------
          !                            ALLOCATE PROBLEM DATA
          !!----------------------------------------------------------------------------------------
           SUBROUTINE allocate_prob_data(nblock, sizeblock)
               !
               ! Allocates the number of blocks 'nblock' and the size of each block 'sizeblock' for the block-separable objective function
               ! 
               ! 
               ! NOTICE: * 'nblock' > 0
               !         * 'sizeblock' > 0               
               !         * 'user_n = nblock * sizeblock
               ! 
               !         * 'user_n' is the dimension of the problem               
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               INTEGER, INTENT(IN) :: nblock             ! the number of blocks in the block-separable objective function
               INTEGER, DIMENSION(nblock), INTENT(IN) :: sizeblock          ! the size of each block in the block-separable objective function           
 
               !**************************** OTHER VARIABLES **************************************
               INTEGER :: i    ! help variable
              
             ! The number of blocks              
               IF (nblock>0) THEN 
                  number_of_blocks = nblock
               ELSE
                  WRITE(*,*) 'There is no blocks in the objective function!'
               END IF 

               ! The size of each block  
               ALLOCATE(size_of_block(number_of_blocks))
               DO i = 1, number_of_blocks
                   IF (sizeblock(i)>0) THEN 
                      size_of_block(i) = sizeblock(i)
                   ELSE
                      WRITE(*,*) 'The size of block ', i, ' cannot be 0 or negative!'
                   END IF 
               END DO
        
           END SUBROUTINE allocate_prob_data
           
          !!----------------------------------------------------------------------------------------
          !                          DEALLOCATE PROBLEM DATA
          !!----------------------------------------------------------------------------------------
           SUBROUTINE deallocate_prob_data()
               !
               ! Deallocates the vector 'sizeblock'
               !               
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************          
 
               !**************************** OTHER VARIABLES **************************************
 
               DEALLOCATE(size_of_block)
        
           END SUBROUTINE deallocate_prob_data
           
        
        !********************************************************************************
        !                                                                               |
        !              FUNCTION VALUES OF THE DC COMPONENTS f_1 AND f_2                 |
        !                                                                               |
        !********************************************************************************
        
           FUNCTION f1(y, problem1, user_n) RESULT(f)       
                !
                ! Calculates the function value of the DC component f_1 at a point 'y'.
                ! Variable 'problem1' identifies the objective function used.
                !
                ! NOTICE: The dimension of 'y' has to be 'user_n'.
                !
                IMPLICIT NONE
                !**************************** NEEDED FROM USER *************************************
                REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: y    ! a point where the function value of the DC component f_1 is calculated
                INTEGER, DIMENSION(number_of_blocks), INTENT(IN) :: problem1                 ! the objective function f_1 for which the value is calculated
                INTEGER, INTENT(IN) :: user_n                   ! the dimension of the problem
                !**************************** OTHER VARIABLES **************************************
                REAL(KIND=dp) :: f                              ! the function value of the DC component f_1 at a point 'y'                                            
                INTEGER :: p                                    ! number_of_blocks
                INTEGER, DIMENSION(number_of_blocks) :: m       ! size_of_block
                INTEGER :: i                                    ! help variables      
                INTEGER :: bstart, bend                         ! help variables      
                
                p = number_of_blocks
                m = size_of_block
                
                f = 0.0_dp
                bstart = 1
                DO i = 1,p
                   bend = bstart+m(i)-1
                   f = f + f1_block(y(bstart:bend), problem1(i), m(i), i)  ! f_1 at 'y'
                   bstart = bstart + m(i)
                END DO
                
           END FUNCTION f1      
           
           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
           
           FUNCTION f2(y, problem2, user_n) RESULT(f)           
                !
                ! Calculates the function value of DC component f_2 at a point 'y'.
                ! Variable 'problem2' identifies the objective function used.
                !
                ! NOTICE: The dimension of 'y' has to be 'user_n'.
                !
                IMPLICIT NONE
                !**************************** NEEDED FROM USER *************************************
                REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: y    ! a point where the function value of the DC component f_2 is calculated
                INTEGER, DIMENSION(number_of_blocks), INTENT(IN) :: problem2                 ! the objective function f_2 for which the value is calculated
                INTEGER, INTENT(IN) :: user_n                   ! the dimension of the problem
                !**************************** OTHER VARIABLES **************************************
                REAL(KIND=dp) :: f                              ! the function value of the DC component f_2 at a point 'y'                                            
                INTEGER :: p                                    ! number_of_blocks
                INTEGER, DIMENSION(number_of_blocks) :: m       ! size_of_block
                INTEGER :: i                                    ! help variables      
                INTEGER :: bstart, bend                         ! help variables      
                
                p = number_of_blocks
                m = size_of_block
                
                f = 0.0_dp
                bstart = 1
                DO i = 1,p
                   bend = bstart+m(i)-1
                   f = f + f2_block(y(bstart:bend), problem2(i), m(i), i)  ! f_2 at 'y'
                   bstart = bstart + m(i)
                END DO

           END FUNCTION f2
           
           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
           
           FUNCTION f1_block(y, problem1, user_n_block, ind) RESULT(f)      
                !
                ! Calculates the function value of the DC component f_1 at a point 'y'.
                ! Variable 'problem1' identifies the objective function used.
                !
                ! NOTICE: The dimension of 'y' has to be 'user_n_block'.
                !
                IMPLICIT NONE
                !**************************** NEEDED FROM USER *************************************
                REAL(KIND=dp), DIMENSION(user_n_block), INTENT(IN) :: y   ! a point where the function value of the DC component f_1 is calculated
                INTEGER, INTENT(IN) :: problem1                           ! the objective function f_1 for which the value is calculated
                INTEGER, INTENT(IN) :: user_n_block                       ! the dimension of the block
                INTEGER, INTENT(IN) :: ind                                ! the index of the block
                !**************************** OTHER VARIABLES **************************************
                REAL(KIND=dp) :: f                              ! the function value of the DC component f_1 at a point 'y'                      
                REAL(KIND=dp) :: apu, apu0, apu1, apu2, apu3    ! help variables                   
                INTEGER :: m                                    ! size_of_block
                INTEGER :: i, l                                 ! help variables  

                REAL(KIND=dp), PARAMETER :: zero=0.0_dp, one=1.0_dp
                REAL(KIND=dp), PARAMETER :: two=2.0_dp              

                m = user_n_block
                
                SELECT CASE(problem1)
                
                   !-------------------------------------
                   ! ************************************ 
                   !           Used in GROUP 0
                   !
                   !      i.e. separable problems
                   ! ************************************ 
                   !-------------------------------------

                   !-------------------------------------
                   !             Function S1
                   !-------------------------------------
                   CASE(1)  
                     f = 0.0_dp
                     
                     DO i = 1, m
                        f = f + y(i)**2
                     END DO
                   !-------------------------------------  
                   
                   !-------------------------------------
                   !             Function S2
                   !-------------------------------------
                   CASE(2)  
                     f = 0.0_dp
                     
                     DO i = 1, m
                        f = f + ABS(y(i))
                     END DO
                   !-------------------------------------    
                  
                   !-------------------------------------
                   !             Function S3
                   !-------------------------------------
                   CASE(3)
                     f = 0.0_dp
                     
                     DO i = 1, m
                        f = f + MAX(0.0_dp, y(i))
                     END DO                 
                   !-------------------------------------                  
                   
                   !-------------------------------------
                   !             Function S4
                   !-------------------------------------
                   CASE(4)
                     f = 0.0_dp
                     
                     DO i = 1, m
                        f = f + MAX(0.0_dp, 1-y(i))
                     END DO                 
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function S5
                   !-------------------------------------
                   CASE(5)
                     f = 0.0_dp
                     
                     DO i = 1, m
                        f = f + MAX(-y(i), EXP(y(i)))
                     END DO 
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function S6
                   !-------------------------------------
                   CASE(6)
                     f = 0.0_dp
                     
                     DO i = 1, m
                        f = f + MAX(y(i), EXP(-y(i)))
                     END DO 
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function S7
                   !-------------------------------------
                   CASE(7)
                     f = 0.0_dp
                     
                     DO i = 1, m
                        f = f + ABS(y(i)) + MAX(0.0_dp, 20.0_dp*(y(i)**2 -y(i) -1.0_dp))
                     END DO 
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function S8
                   !-------------------------------------
                   CASE(8)
                     f = 0.0_dp
                     
                     DO i = 1, m
                        f = f + MAX(y(i)+1.0_dp, -2.0_dp*y(i)-2.0_dp, 0.5_dp*y(i)+5.0_dp)
                     END DO
                     
                   !-------------------------------------
                   
                   !-------------------------------------
                   ! ************************************ 
                   !           Used in GROUP 1
                   !
                   !   i.e. block-separable problems
                   ! ************************************ 
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function B1
                   !-------------------------------------
                   CASE(1001)             
                     
                     apu = -100.0_dp
                     DO i = 1, m
                        IF (apu < y(i)**2 ) THEN 
                            apu = y(i)**2
                        END IF
                     END DO
                     f = (m+1)*apu
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function B2
                   !-------------------------------------
                   CASE(1002)             
                     
                     apu = -100.0_dp
                     DO i = 1, m
                        IF (apu < y(i)**2 ) THEN 
                            apu = y(i)**2
                        END IF
                     END DO
                     f = apu
                   !-------------------------------------

                   !-------------------------------------
                   !             Function B3
                   !-------------------------------------
                   CASE(1003)             
                     f = 0.0_dp
                     
                     DO i = 1, m
                        f = f + y(i)**2
                     END DO
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function B4
                   !-------------------------------------
                   CASE(1004)             
                     
                     apu = -100.0_dp
                     DO i = 1, m
                        IF ( apu < ABS(y(i)) ) THEN 
                            apu = ABS(y(i))
                        END IF
                     END DO
                     f = m*apu
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function B5
                   !-------------------------------------
                   CASE(1005)             
                     f = 0.0_dp
                     
                     DO i = 1, m
                        f = f + ABS(y(i))
                     END DO
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function B6
                   !-------------------------------------
                   CASE(1006)             
                     f = 0.0_dp
                     
                     DO i = 2, m
                        f = f + ABS(y(i)-y(i-1))
                     END DO
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function B7
                   !-------------------------------------
                   CASE(1007)             
                     
                     apu1 = -100.0_dp
                     DO l = 1, m
                        apu2 = 0.0_dp
                        DO i = 1, m
                            apu2 = apu2 + (y(i) / REAL(l+i-1, KIND=dp))
                        END DO
                        IF (apu1 < ABS(apu2)) THEN
                            apu1 = ABS(apu2)
                        END IF
                     END DO
                     f = apu1
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function B8
                   !-------------------------------------
                   CASE(1008)             
                     f = 0.0_dp
                     
                     DO i = 1, m-1
                        apu1 = -y(i)-y(i+1)
                        apu2 = -y(i)-y(i+1)+y(i)*y(i)+y(i+1)*y(i+1)-1.0_dp
                        IF (apu1 >= apu2) THEN
                            f = f + apu1
                        ELSE
                            f = f + apu2
                        END IF
                    END DO
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function B9
                   !-------------------------------------
                   CASE(1009)             
                     f = 0.0_dp
                     
                     DO i = 1, m-1
                        apu1 = y(i)*y(i)*y(i)*y(i)+y(i+1)*y(i+1)
                        apu2 = (2.0_dp-y(i))*(2.0_dp-y(i))+(2.0_dp-y(i+1))*(2.0_dp-y(i+1))
                        apu3 = 2.0_dp*EXP(-y(i)+y(i+1))
                        apu = MAX(apu1, apu2, apu3)
                        f = f + apu
                     END DO
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function B10
                   !-------------------------------------
                   CASE(1010)             

                     apu1 = 0.0_dp
                     apu2 = 0.0_dp
                     apu3 = 0.0_dp
                     DO i = 1, m-1
                        apu1 = apu1 + y(i)*y(i)*y(i)*y(i)+y(i+1)*y(i+1)
                        apu2 = apu2 + (2.0_dp-y(i))*(2.0_dp-y(i))+(2.0_dp-y(i+1))*(2.0_dp-y(i+1))
                        apu3 = apu3 + 2.0_dp*EXP(-y(i)+y(i+1))
                     END DO
                     f = MAX(apu1, apu2, apu3)
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function B11
                   !-------------------------------------
                   CASE(1011)             
                     f = 0.0_dp
                     
                     DO i = 1, m
                        f = f + ABS(y(i)) 
                        f = f + MAX(20.0_dp*(y(i)*y(i)-y(i)-1.0_dp), 0.0_dp)
                     END DO
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function B12
                   !-------------------------------------
                   CASE(1012)             
                     
                     apu1 = -100.0_dp
                     DO l = 1, m
                        apu2 = 0.0_dp
                        DO i = 1, m
                            apu2 = apu2 + (y(i) / REAL(l+i-1, KIND=dp))
                        END DO
                        IF (apu1 < ABS(apu2)) THEN
                            apu1 = ABS(apu2)
                        END IF
                     END DO
                     f = m*apu1
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function B13
                   !-------------------------------------
                   CASE(1013)             
                     f = 0.0_dp
                     
                     DO l = 1, m
                        apu = 0.0_dp
                        DO i = 1, m
                            apu = apu + (y(i) / REAL(l+i-1, KIND=dp))
                        END DO
                        f = f + ABS(apu)
                     END DO
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function B14
                   !-------------------------------------
                   CASE(1014)             
                     
                     apu0 = -100.0_dp
                     DO i = 1, m-1
                        apu1 = y(i)*y(i)*y(i)*y(i)+y(i+1)*y(i+1)
                        apu2 = (2.0_dp-y(i))*(2.0_dp-y(i))+(2.0_dp-y(i+1))*(2.0_dp-y(i+1))
                        apu3 = 2.0_dp*EXP(-y(i)+y(i+1))
                        apu = MAX(apu1, apu2, apu3)
                        IF (apu0 < apu) THEN
                            apu0 = apu
                        END IF
                     END DO
                     f = (m-1)*apu0
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function B15
                   !-------------------------------------
                   CASE(1015)             
                     
                     apu = 0.0_dp
                     DO i = 1, m-1
                        apu = apu + y(i)**2 + (y(i+1)-1.0_dp)**2
                        apu = apu + y(i+1) - 1.0_dp
                     END DO
                     apu = 2.0_dp * apu
                     IF (apu > 0.0_dp) THEN
                        f = apu
                     END IF 
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function B16
                   !-------------------------------------
                   CASE(1016)             
                     f = 0.0_dp
                     
                     DO i = 1, m-1
                        f = f + y(i)**2 + (y(i+1)-1.0_dp)**2
                        f = f + y(i+1) - 1.0_dp
                     END DO
                   !-------------------------------------
                            
                END SELECT
                
           END FUNCTION f1_block      
           
           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                   
           FUNCTION f2_block(y, problem2, user_n_block, ind) RESULT(f)           
                !
                ! Calculates the function value of DC component f_2 at a point 'y'.
                ! Variable 'problem2' identifies the objective function used.
                !
                ! NOTICE: The dimension of 'y' has to be 'user_n_block'.
                !
                IMPLICIT NONE
                !**************************** NEEDED FROM USER *************************************
                REAL(KIND=dp), DIMENSION(user_n_block), INTENT(IN) :: y   ! a point where the function value of the DC component f_2 is calculated
                INTEGER, INTENT(IN) :: problem2                           ! the objective function f_2 for which the value is calculated   
                INTEGER, INTENT(IN) :: user_n_block                       ! the dimension of the block             
                INTEGER, INTENT(IN) :: ind                                ! the index of the block             
                !**************************** OTHER VARIABLES **************************************
                REAL(KIND=dp) :: f                            ! the function value of the DC component f_2 at a point 'y'
                REAL(KIND=dp) :: apu, apu0, apu1, apu2, apu3  ! help variables 
                INTEGER :: m                                  ! size_of_block
                INTEGER :: i, l                               ! help variables
                
                REAL(KIND=dp), PARAMETER :: zero=0.0_dp, one=1.0_dp
                REAL(KIND=dp), PARAMETER :: two=2.0_dp                 
                
                
                m = user_n_block
                
                SELECT CASE(problem2)
                   
                   !-------------------------------------
                   ! ************************************ 
                   !           Used in GROUP 0
                   !
                   !      i.e. separable problems
                   ! ************************************ 
                   !-------------------------------------

                   !-------------------------------------
                   !             Function S1
                   !-------------------------------------
                   CASE(1)  
                     f = 0.0_dp
                     
                     DO i = 1, m
                        f = f + y(i)**2
                     END DO
                   !-------------------------------------  
                   
                   !-------------------------------------
                   !             Function S2
                   !-------------------------------------
                   CASE(2)  
                     f = 0.0_dp
                     
                     DO i = 1, m
                        f = f + ABS(y(i))
                     END DO
                   !-------------------------------------    
                  
                   !-------------------------------------
                   !             Function S3
                   !-------------------------------------
                   CASE(3)
                     f = 0.0_dp
                     
                     DO i = 1, m
                        f = f + MAX(0.0_dp, y(i))
                     END DO                 
                   !-------------------------------------                  
                   
                   !-------------------------------------
                   !             Function S4
                   !-------------------------------------
                   CASE(4)
                     f = 0.0_dp
                     
                     DO i = 1, m
                        f = f + MAX(0.0_dp, 1-y(i))
                     END DO                 
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function S5
                   !-------------------------------------
                   CASE(5)
                     f = 0.0_dp
                     
                     DO i = 1, m
                        f = f + MAX(-y(i), EXP(y(i)))
                     END DO 
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function S6
                   !-------------------------------------
                   CASE(6)
                     f = 0.0_dp
                     
                     DO i = 1, m
                        f = f + MAX(y(i), EXP(-y(i)))
                     END DO 
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function S7
                   !-------------------------------------
                   CASE(7)
                     f = 0.0_dp
                     
                     DO i = 1, m
                        f = f + ABS(y(i)) + MAX(0.0_dp, 20.0_dp*(y(i)**2 -y(i) -1.0_dp))
                     END DO 
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function S8
                   !-------------------------------------
                   CASE(8)
                     f = 0.0_dp
                     
                     DO i = 1, m
                        f = f + MAX(y(i)+1.0_dp, -2.0_dp*y(i)-2.0_dp, 0.5_dp*y(i)+5.0_dp)
                     END DO
                     
                   !-------------------------------------
                   
                   !-------------------------------------
                   ! ************************************ 
                   !           Used in GROUP 1
                   !
                   !   i.e. block-separable problems
                   ! ************************************ 
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function B1
                   !-------------------------------------
                   CASE(1001)             
                     
                     apu = -100.0_dp
                     DO i = 1, m
                        IF (apu < y(i)**2 ) THEN 
                            apu = y(i)**2
                        END IF
                     END DO
                     f = (m+1)*apu
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function B2
                   !-------------------------------------
                   CASE(1002)             
                     
                     apu = -100.0_dp
                     DO i = 1, m
                        IF (apu < y(i)**2 ) THEN 
                            apu = y(i)**2
                        END IF
                     END DO
                     f = apu
                   !-------------------------------------

                   !-------------------------------------
                   !             Function B3
                   !-------------------------------------
                   CASE(1003)             
                     f = 0.0_dp
                     
                     DO i = 1, m
                        f = f + y(i)**2
                     END DO
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function B4
                   !-------------------------------------
                   CASE(1004)             
                     
                     apu = -100.0_dp
                     DO i = 1, m
                        IF ( apu < ABS(y(i)) ) THEN 
                            apu = ABS(y(i))
                        END IF
                     END DO
                     f = m*apu
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function B5
                   !-------------------------------------
                   CASE(1005)             
                     f = 0.0_dp
                     
                     DO i = 1, m
                        f = f + ABS(y(i))
                     END DO
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function B6
                   !-------------------------------------
                   CASE(1006)             
                     f = 0.0_dp
                     
                     DO i = 2, m
                        f = f + ABS(y(i)-y(i-1))
                     END DO
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function B7
                   !-------------------------------------
                   CASE(1007)             
                     
                     apu1 = -100.0_dp
                     DO l = 1, m
                        apu2 = 0.0_dp
                        DO i = 1, m
                            apu2 = apu2 + (y(i) / REAL(l+i-1, KIND=dp))
                        END DO
                        IF (apu1 < ABS(apu2)) THEN
                            apu1 = ABS(apu2)
                        END IF
                     END DO
                     f = apu1
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function B8
                   !-------------------------------------
                   CASE(1008)             
                     f = 0.0_dp
                     
                     DO i = 1, m-1
                        apu1 = -y(i)-y(i+1)
                        apu2 = -y(i)-y(i+1)+y(i)*y(i)+y(i+1)*y(i+1)-1.0_dp
                        IF (apu1 >= apu2) THEN
                            f = f + apu1
                        ELSE
                            f = f + apu2
                        END IF
                    END DO
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function B9
                   !-------------------------------------
                   CASE(1009)             
                     f = 0.0_dp
                     
                     DO i = 1, m-1
                        apu1 = y(i)*y(i)*y(i)*y(i)+y(i+1)*y(i+1)
                        apu2 = (2.0_dp-y(i))*(2.0_dp-y(i))+(2.0_dp-y(i+1))*(2.0_dp-y(i+1))
                        apu3 = 2.0_dp*EXP(-y(i)+y(i+1))
                        apu = MAX(apu1, apu2, apu3)
                        f = f + apu
                     END DO
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function B10
                   !-------------------------------------
                   CASE(1010)             

                     apu1 = 0.0_dp
                     apu2 = 0.0_dp
                     apu3 = 0.0_dp
                     DO i = 1, m-1
                        apu1 = apu1 + y(i)*y(i)*y(i)*y(i)+y(i+1)*y(i+1)
                        apu2 = apu2 + (2.0_dp-y(i))*(2.0_dp-y(i))+(2.0_dp-y(i+1))*(2.0_dp-y(i+1))
                        apu3 = apu3 + 2.0_dp*EXP(-y(i)+y(i+1))
                     END DO
                     f = MAX(apu1, apu2, apu3)
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function B11
                   !-------------------------------------
                   CASE(1011)             
                     f = 0.0_dp
                     
                     DO i = 1, m
                        f = f + ABS(y(i)) 
                        f = f + MAX(20.0_dp*(y(i)*y(i)-y(i)-1.0_dp), 0.0_dp)
                     END DO
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function B12
                   !-------------------------------------
                   CASE(1012)             
                     
                     apu1 = -100.0_dp
                     DO l = 1, m
                        apu2 = 0.0_dp
                        DO i = 1, m
                            apu2 = apu2 + (y(i) / REAL(l+i-1, KIND=dp))
                        END DO
                        IF (apu1 < ABS(apu2)) THEN
                            apu1 = ABS(apu2)
                        END IF
                     END DO
                     f = m*apu1
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function B13
                   !-------------------------------------
                   CASE(1013)             
                     f = 0.0_dp
                     
                     DO l = 1, m
                        apu = 0.0_dp
                        DO i = 1, m
                            apu = apu + (y(i) / REAL(l+i-1, KIND=dp))
                        END DO
                        f = f + ABS(apu)
                     END DO
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function B14
                   !-------------------------------------
                   CASE(1014)             
                     
                     apu0 = -100.0_dp
                     DO i = 1, m-1
                        apu1 = y(i)*y(i)*y(i)*y(i)+y(i+1)*y(i+1)
                        apu2 = (2.0_dp-y(i))*(2.0_dp-y(i))+(2.0_dp-y(i+1))*(2.0_dp-y(i+1))
                        apu3 = 2.0_dp*EXP(-y(i)+y(i+1))
                        apu = MAX(apu1, apu2, apu3)
                        IF (apu0 < apu) THEN
                            apu0 = apu
                        END IF
                     END DO
                     f = (m-1)*apu0
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function B15
                   !-------------------------------------
                   CASE(1015)             
                     
                     apu = 0.0_dp
                     DO i = 1, m-1
                        apu = apu + y(i)**2 + (y(i+1)-1.0_dp)**2
                        apu = apu + y(i+1) - 1.0_dp
                     END DO
                     apu = 2.0_dp * apu
                     IF (apu > 0.0_dp) THEN
                        f = apu
                     END IF 
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function B16
                   !-------------------------------------
                   CASE(1016)             
                     f = 0.0_dp
                     
                     DO i = 1, m-1
                        f = f + y(i)**2 + (y(i+1)-1.0_dp)**2
                        f = f + y(i+1) - 1.0_dp
                     END DO
                   !-------------------------------------
                            
                END SELECT
           
           
           END FUNCTION f2_block
           
           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -           


        !********************************************************************************
        !                                                                               |
        !                SUBGRADIENTS OF THE DC COMPONENTS f_1 AND f_2                  |
        !                                                                               |       
        !********************************************************************************       
        
           FUNCTION subgradient_f1(y, problem1, user_n) RESULT(grad)
                !
                ! Calculates a subgradient of the DC component f_1 at a point 'y'.
                ! Variable 'problem1' identifies the objective function used.
                !
                ! NOTICE: * The dimension of 'y' has to be 'user_n'.
                !         * The dimension of 'grad' is 'user_n'.
                !
                IMPLICIT NONE
                !**************************** NEEDED FROM USER *************************************
                REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: y    ! a point where the subgradient of the DC component f_1 is calculated
                INTEGER, DIMENSION(number_of_blocks), INTENT(IN) :: problem1                 ! the objective function f_1 for which the subgradient is calculated      
                INTEGER, INTENT(IN) :: user_n                   ! the dimension of the problem              
                !**************************** OTHER VARIABLES **************************************
                REAL(KIND=dp), DIMENSION(SIZE(y)) :: grad       ! the subgradient of the DC component f_1 at a point 'y'            
                INTEGER :: p                                    ! number_of_blocks
                INTEGER, DIMENSION(number_of_blocks) :: m       ! size_of_block
                INTEGER :: bstart                               ! the first element of block
                INTEGER :: bend                                 ! the last element of block
                INTEGER :: i
                
                p = number_of_blocks
                m = size_of_block
                
                bstart = 1
                DO i = 1,p
                    bend = bstart+m(i)-1
                    grad(bstart:bend) = subgradient_f1_block(y(bstart:bend), problem1(i), m(i), i)  ! Subgradient of f_1 at 'y'
                    bstart = bstart + m(i)
                END DO                
                
           END FUNCTION subgradient_f1      
           
           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
           
           FUNCTION subgradient_f2(y, problem2, user_n) RESULT(grad)                
                !
                ! Calculate a subgradient of the DC component f_2 at a point 'y'.
                ! Variable 'problem2' identifies the objective function used.
                !
                ! NOTICE: * The dimension of 'y' has to be 'user_n'.
                !         * The dimension of 'grad' is also 'user_n'.
                !
                IMPLICIT NONE
                !**************************** NEEDED FROM USER *************************************
                REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: y    ! a point where the subgradient of the DC component f_2 is calculated
                INTEGER, DIMENSION(number_of_blocks), INTENT(IN) :: problem2                 ! the objective function f_2 for which the subgradient is calculated      
                INTEGER, INTENT(IN) :: user_n                   ! the dimension of the problem              
                !**************************** OTHER VARIABLES **************************************
                REAL(KIND=dp), DIMENSION(SIZE(y)) :: grad       ! the subgradient of the DC component f_2 at a point 'y'            
                INTEGER :: p                                    ! number_of_blocks
                INTEGER, DIMENSION(number_of_blocks) :: m       ! size_of_block
                INTEGER :: bstart                               ! the first element of block
                INTEGER :: bend                                 ! the last element of block
                INTEGER :: i
                
                p = number_of_blocks
                m = size_of_block
                
                bstart = 1
                DO i = 1,p
                    bend = bstart+m(i)-1
                    grad(bstart:bend) = subgradient_f2_block(y(bstart:bend), problem2(i), m(i), i)  ! Subgradient of f_2 at 'y'
                    bstart = bstart + m(i)
                END DO    

           END FUNCTION subgradient_f2

           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -      
           
           FUNCTION subgradient_f1_block(y, problem1, user_n_block, indb) RESULT(grad)
                !
                ! Calculates a subgradient of the DC component f_1 at a point 'y'.
                ! Variable 'problem1' identifies the objective function used.
                !
                ! NOTICE: * The dimension of 'y' has to be 'user_n_block'.
                !         * The dimension of 'grad' is 'user_n_block'.
                !
                IMPLICIT NONE
                !**************************** NEEDED FROM USER *************************************
                REAL(KIND=dp), DIMENSION(user_n_block), INTENT(IN) :: y   ! a point where the subgradient of the DC component f_1 is calculated
                INTEGER, INTENT(IN) :: problem1                           ! the objective function f_1 for which the subgradient is calculated      
                INTEGER, INTENT(IN) :: user_n_block                       ! the dimension of the block     
                INTEGER, INTENT(IN) :: indb                               ! the index of the block              
                !**************************** OTHER VARIABLES **************************************
                REAL(KIND=dp), DIMENSION(SIZE(y)) :: grad       ! the subgradient of the DC component f_1 at a point 'y'           
                REAL(KIND=dp) :: apu, apu0, apu1, apu2, apu3    ! help variables
                INTEGER :: m                                    ! size_of_block
                INTEGER :: i, l, ind, ind2                      ! help variables

                REAL(KIND=dp), PARAMETER :: zero=0.0_dp, one=1.0_dp
                REAL(KIND=dp), PARAMETER :: two=2.0_dp              
                
                REAL(KIND=dp) :: f_sign1                        ! sign of the active partial objective for f1 and for f2                                  

                
                m = user_n_block

                SELECT CASE(problem1)
                   
                   !-------------------------------------
                   ! ************************************ 
                   !           Used in GROUP 0
                   !
                   !      i.e. separable problems
                   ! ************************************ 
                   !-------------------------------------

                   !-------------------------------------
                   !             Function S1
                   !-------------------------------------
                   CASE(1)  
                     
                     grad = 2.0_dp*y
                   !-------------------------------------  
                   
                   !-------------------------------------
                   !             Function S2
                   !-------------------------------------
                   CASE(2)  
                     
                     DO i = 1, m
                        IF (y(i) <= 0.0_dp) THEN
                            grad(i) = -1.0_dp
                        ELSE
                            grad(i) = 1.0_dp
                        END IF
                     END DO
                   !-------------------------------------                  
                  
                   !-------------------------------------
                   !             Function S3
                   !-------------------------------------
                   CASE(3)
                     grad = 0.0_dp
                     
                     DO i = 1, m
                        IF (y(i) > 0.0_dp) THEN
                            grad(i) = 1.0_dp
                        END IF
                     END DO                 
                   !-------------------------------------                      
                   
                   !-------------------------------------
                   !             Function S4
                   !-------------------------------------
                   CASE(4)
                     grad = 0.0_dp
                     
                     DO i = 1, m
                        IF (1.0_dp-y(i) > 0.0_dp) THEN
                            grad(i) = -1.0_dp
                        END IF
                     END DO                 
                   !-------------------------------------   
                   
                   !-------------------------------------
                   !             Function S5
                   !-------------------------------------
                   CASE(5)
                     
                     DO i = 1, m
                        IF (EXP(y(i)) >= -y(i)) THEN
                            grad(i) = EXP(y(i))
                        ELSE
                            grad(i) = -1.0_dp
                        END IF
                     END DO 
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function S6
                   !-------------------------------------
                   CASE(6)
                     
                     DO i = 1, m
                        IF (EXP(-y(i)) > y(i)) THEN
                            grad(i) = -EXP(-y(i))
                        ELSE
                            grad(i) = 1.0_dp
                        END IF
                     END DO 
                   !-------------------------------------    
                   
                   !-------------------------------------
                   !             Function S7
                   !-------------------------------------
                   CASE(7)
                     
                     DO i = 1, m
                        IF (y(i) > 0.0_dp) THEN
                            grad(i) = 1.0_dp
                        ELSE    
                            grad(i) = -1.0_dp
                        END IF
                        IF (20.0_dp*(y(i)**2-y(i)-1.0_dp) > 0.0_dp) THEN
                            grad(i) = grad(i) + 40.0_dp*y(i)-20.0_dp
                        END IF
                     END DO 
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function S8
                   !-------------------------------------
                   CASE(8)
                     
                     DO i = 1, m
                        apu1 = y(i)+1.0_dp
                        apu2 = -2.0_dp*y(i)-2.0_dp
                        apu3 = 0.5_dp*y(i)+5.0_dp
                        apu = MAX(apu1, apu2, apu3)
                        IF (apu == apu1) THEN
                            grad(i) = 1.0_dp
                        ELSE IF (apu == apu2) THEN
                            grad(i) = -2.0_dp
                        ELSE
                            grad(i) = 0.5_dp
                        END IF
                     END DO
                   !-------------------------------------
                   
                   !-------------------------------------
                   ! ************************************ 
                   !          Used in GROUP 1
                   !
                   !   i.e. block-separable problems
                   ! ************************************ 
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function B1
                   !-------------------------------------
                   CASE(1001)             
                     grad = 0.0_dp
                     ind = 1
                     
                     apu = -100.0_dp
                     DO i = 1, m
                        IF (apu < y(i)**2 ) THEN 
                            apu = y(i)**2
                            ind = i
                        END IF
                     END DO
                     grad(ind) = grad(ind) + 2.0_dp*(m+1)*y(ind)
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function B2
                   !-------------------------------------
                   CASE(1002) 
                     grad = 0.0_dp
                     ind = 1
                     
                     apu = -100.0_dp
                     DO i = 1, m
                        IF (apu < y(i)**2) THEN 
                            apu = y(i)**2
                            ind = i
                        END IF
                     END DO
                     grad(ind) = grad(ind) + 2.0_dp*y(ind)
                   !-------------------------------------              

                   !-------------------------------------
                   !             Function B3
                   !-------------------------------------
                   CASE(1003)    
                     grad = 2.0_dp * y
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function B4
                   !-------------------------------------
                   CASE(1004) 
                     grad = 0.0_dp
                     ind = 1
                     
                     apu = -100.0_dp
                     DO i = 1, m
                        IF (apu < ABS(y(i))) THEN 
                            apu = ABS(y(i))
                            ind = i
                        END IF
                     END DO
                     IF (y(ind) <= 0.0_dp) THEN 
                        grad(ind) = grad(ind) - m
                     ELSE 
                        grad(ind) = grad(ind) + m
                     END IF 
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function B5
                   !-------------------------------------
                   CASE(1005)             
                     grad = 0.0_dp
                     
                     DO i = 1, m
                        IF (y(i) <= 0.0_dp) THEN 
                            grad(i) = grad(i) - 1.0_dp
                        ELSE
                            grad(i) = grad(i) + 1.0_dp
                        END IF
                     END DO
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function B6
                   !-------------------------------------
                   CASE(1006)             
                     grad = 0.0_dp
                     
                     DO i = 2, m
                        IF ( y(i) - y(i-1) <= 0.0_dp) THEN
                            grad(i-1) = grad(i-1) + 1.0_dp 
                            grad(i) = grad(i) - 1.0_dp 
                        ELSE
                            grad(i-1) = grad(i-1) - 1.0_dp  
                            grad(i) = grad(i) + 1.0_dp     
                        END IF
                     END DO
                   !-------------------------------------   
                   
                   !-------------------------------------
                   !             Function B7
                   !-------------------------------------
                   CASE(1007)
                     grad = 0.0_dp
                     
                     apu1 = -100.0_dp
                     DO l = 1, m
                        apu2 = 0.0_dp
                        DO i = 1, m
                            apu2 = apu2 + (y(i) / REAL(l+i-1, KIND=dp))
                        END DO
                        IF (apu1 < ABS(apu2)) THEN
                            apu1 = ABS(apu2)
                            ind = l
                            f_sign1 = SIGN(1.0_dp,apu2)
                        END IF
                     END DO
                     DO i = 1, m
                        grad(i) = grad(i)+f_sign1/REAL(ind+i-1, KIND=dp)
                     END DO 
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function B8
                   !-------------------------------------
                   CASE(1008)             
                     grad = 0.0_dp
                     
                     DO i = 1, m-1
                        apu1 = -y(i)-y(i+1)
                        apu2 = -y(i)-y(i+1)+y(i)*y(i)+y(i+1)*y(i+1)-1.0_dp
                        IF (apu1 >= apu2) THEN
                            grad(i) = grad(i) - 1.0_dp
                            grad(i+1) = grad(i+1) - 1.0_dp
                        ELSE
                            grad(i) = grad(i)-1.0_dp+2.0_dp*y(i)
                            grad(i+1) = grad(i+1)-1.0_dp+2.0_dp*y(i+1)
                        END IF
                     END DO
                   !-------------------------------------

                   !-------------------------------------
                   !             Function B9
                   !-------------------------------------
                   CASE(1009)             
                     grad = 0.0_dp
                     
                     DO i = 1, m-1
                        apu1 = y(i)*y(i)*y(i)*y(i)+y(i+1)*y(i+1)
                        apu2 = (2.0_dp-y(i))*(2.0_dp-y(i))+(2.0_dp-y(i+1))*(2.0_dp-y(i+1))
                        apu3 = 2.0_dp*EXP(-y(i)+y(i+1))
                        apu = MAX(apu1, apu2, apu3)
                        IF (apu == apu1) THEN
                            grad(i) = grad(i)+4.0_dp*y(i)*y(i)*y(i)
                            grad(i+1) = grad(i+1)+2.0_dp*y(i+1)
                        ELSE IF (apu == apu2) THEN
                            grad(i) = grad(i)+2.0_dp*y(i)-4.0_dp
                            grad(i+1) = grad(i+1)+2.0_dp*y(i+1)-4.0_dp
                        ELSE
                            grad(i) = grad(i) - 2.0_dp*EXP(-y(i)+y(i+1))
                            grad(i+1) = grad(i+1) + 2.0_dp*EXP(-y(i)+y(i+1))
                        END IF   
                     END DO
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function B10
                   !-------------------------------------
                   CASE(1010)             
                     grad = 0.0_dp
                     
                     apu1 = 0.0_dp
                     apu2 = 0.0_dp
                     apu3 = 0.0_dp
                     DO i = 1, m-1
                        apu1 = apu1 + y(i)*y(i)*y(i)*y(i)+y(i+1)*y(i+1)
                        apu2 = apu2 + (2.0_dp-y(i))*(2.0_dp-y(i))+(2.0_dp-y(i+1))*(2.0_dp-y(i+1))
                        apu3 = apu3 + 2.0_dp*EXP(-y(i)+y(i+1))
                     END DO
                     apu = MAX(apu1, apu2, apu3)
                     IF (apu == apu1) THEN
                        DO i = 1, m-1
                            grad(i) = grad(i)+4.0_dp*y(i)*y(i)*y(i)
                            grad(i+1) = grad(i+1)+2.0_dp*y(i+1)
                        END DO 
                     ELSE IF (apu == apu2) THEN
                        DO i = 1, m-1
                            grad(i) = grad(i)+2.0_dp*y(i)-4.0_dp
                            grad(i+1) = grad(i+1)+2.0_dp*y(i+1)-4.0_dp
                        END DO
                     ELSE
                        DO i = 1, m-1
                            grad(i) = grad(i) - 2.0_dp*EXP(-y(i)+y(i+1))
                            grad(i+1) = grad(i+1) + 2.0_dp*EXP(-y(i)+y(i+1))
                        END DO
                     END IF
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function B11
                   !-------------------------------------
                   CASE(1011)             
                     grad = 0.0_dp
                     
                     DO i = 1, m
                        IF ( y(i) < 0.0_dp ) THEN 
                           grad(i) = grad(i)-1.0_dp
                        ELSE
                           grad(i) = grad(i)+1.0_dp
                        END IF
                        IF (20.0_dp*(y(i)*y(i)-y(i)-1.0_dp) > 0.0_dp) THEN
                           grad(i) = grad(i) + 20.0_dp*(2.0_dp*y(i)-1.0_dp)
                        END IF
                     END DO
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function B12
                   !-------------------------------------
                   CASE(1012)
                     grad = 0.0_dp
                     
                     apu1 = -100.0_dp
                     DO l = 1, m
                        apu2 = 0.0_dp
                        DO i = 1, m
                            apu2 = apu2 + (y(i) / REAL(l+i-1, KIND=dp))
                        END DO
                        IF (apu1 < ABS(apu2)) THEN
                            apu1 = ABS(apu2)
                            ind = l
                            f_sign1 = SIGN(1.0_dp,apu2)
                        END IF
                     END DO
                     DO i = 1, m
                        grad(i) = grad(i)+m*f_sign1/REAL(ind+i-1, KIND=dp)
                     END DO 
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function B13
                   !-------------------------------------
                   CASE(1013)             
                     grad = 0.0_dp
                     
                     DO l = 1, m
                        apu = 0.0_dp
                        DO i = 1, m
                            apu = apu + y(i) / REAL(l+i-1, KIND=dp)
                        END DO
                        f_sign1 = SIGN(1.0_dp,apu)
                        DO i = 1, m
                            grad(i) = grad(i) + f_sign1/REAL(l+i-1, KIND=dp)
                        END DO
                     END DO
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function B14
                   !-------------------------------------
                   CASE(1014)             
                     grad = 0.0_dp

                     apu0 = -100.0_dp
                     DO i = 1, m-1
                        apu1 = y(i)*y(i)*y(i)*y(i)+y(i+1)*y(i+1)
                        apu2 = (2.0_dp-y(i))*(2.0_dp-y(i))+(2.0_dp-y(i+1))*(2.0_dp-y(i+1))
                        apu3 = 2.0_dp*EXP(-y(i)+y(i+1))
                        apu = MAX(apu1, apu2, apu3)
                        IF (apu0 < apu) THEN
                            apu0 = apu
                            ind = i
                            IF (apu == apu1) THEN
                                ind2 = 1
                            ELSE IF (apu == apu2) THEN
                                ind2 = 2
                            ELSE
                                ind2 = 3
                            END IF
                        END IF
                     END DO
                     IF (ind2 == 1) THEN
                        grad(ind) = grad(ind)+(m-1)*(4.0_dp*y(ind)*y(ind)*y(ind))
                        grad(ind+1) = grad(ind+1)+(m-1)*2.0_dp*y(ind+1)
                     ELSE IF (ind2 == 2) THEN
                        grad(ind) = grad(ind)+(m-1)*(2.0_dp*y(ind)-4.0_dp)
                        grad(ind+1) = grad(ind+1)+(m-1)*(2.0_dp*y(ind+1)-4.0_dp)
                     ELSE
                        grad(ind) = grad(ind) - (m-1)*(2.0_dp*EXP(-y(ind)+y(ind+1)))
                        grad(ind+1) = grad(ind+1) + (m-1)*(2.0_dp*EXP(-y(ind)+y(ind+1)))
                     END IF
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function B15
                   !-------------------------------------
                   CASE(1015)             
                     grad = 0.0_dp
                     
                     apu = 0.0_dp
                     DO i = 1, m-1
                        apu = apu + y(i)**2 + (y(i+1)-1.0_dp)**2
                        apu = apu + y(i+1) - 1.0_dp
                     END DO
                     apu = 2.0_dp * apu
                     IF (apu > 0.0_dp) THEN
                        DO i = 1, m-1
                            grad(i) = grad(i) + 4.0_dp *y(i)
                            grad(i+1) = grad(i+1) + 4.0_dp *(y(i+1)-1.0_dp) +2.0_dp
                        END DO
                     END IF 
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function B16
                   !-------------------------------------
                   CASE(1016)             
                     grad = 0.0_dp
                     
                     DO i = 1, m-1
                        grad(i) = grad(i) + 2.0_dp * y(i) 
                        grad(i+1) = grad(i+1) + 2.0_dp * (y(i+1) -1.0_dp) +1.0_dp
                     END DO
                   !-------------------------------------

                END SELECT                      
                
           END FUNCTION subgradient_f1_block      
           
           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

           FUNCTION subgradient_f2_block(y, problem2, user_n_block, indb) RESULT(grad)                
                !
                ! Calculate a subgradient of the DC component f_2 at a point 'y'.
                ! Variable 'problem2' identifies the objective function used.
                !
                ! NOTICE: * The dimension of 'y' has to be 'user_n_block'.
                !         * The dimension of 'grad' is also 'user_n_block'.
                !
                IMPLICIT NONE
                !**************************** NEEDED FROM USER *************************************
                REAL(KIND=dp), DIMENSION(user_n_block), INTENT(IN) :: y   ! a point where the subgradient of the DC component f_2 is calculated
                INTEGER, INTENT(IN) :: problem2                           ! the objective function f_2 for which the subgradient is calculated
                INTEGER, INTENT(IN) :: user_n_block                       ! the dimension of the block 
                INTEGER, INTENT(IN) :: indb                               ! the index of the block              
                !**************************** OTHER VARIABLES **************************************
                REAL(KIND=dp), DIMENSION(SIZE(y)) :: grad       ! the subgradient of the DC component f_2 at a point 'y'    
                REAL(KIND=dp) :: apu, apu0, apu1, apu2, apu3    ! help varaibles
                INTEGER :: m                                    ! size_of_block
                INTEGER :: i, l, ind, ind2                      ! help variables

                REAL(KIND=dp), PARAMETER :: zero=0.0_dp, one=1.0_dp
                REAL(KIND=dp), PARAMETER :: two=2.0_dp              
                
                REAL(KIND=dp) :: f_sign1                        ! sign of the active partial objective for f1 and for f2                                     
               
                
                m = user_n_block
                
                SELECT CASE(problem2)           
                   
                   !-------------------------------------
                   ! ************************************ 
                   !           Used in GROUP 0
                   !
                   !      i.e. separable problems
                   ! ************************************ 
                   !-------------------------------------

                   !-------------------------------------
                   !             Function S1
                   !-------------------------------------
                   CASE(1)  
                     
                     grad = 2.0_dp*y
                   !-------------------------------------  
                   
                   !-------------------------------------
                   !             Function S2
                   !-------------------------------------
                   CASE(2)  
                     
                     DO i = 1, m
                        IF (y(i) <= 0.0_dp) THEN
                            grad(i) = -1.0_dp
                        ELSE
                            grad(i) = 1.0_dp
                        END IF
                     END DO
                   !-------------------------------------                  
                  
                   !-------------------------------------
                   !             Function S3
                   !-------------------------------------
                   CASE(3)
                     grad = 0.0_dp
                     
                     DO i = 1, m
                        IF (y(i) > 0.0_dp) THEN
                            grad(i) = 1.0_dp
                        END IF
                     END DO                 
                   !-------------------------------------                      
                   
                   !-------------------------------------
                   !             Function S4
                   !-------------------------------------
                   CASE(4)
                     grad = 0.0_dp
                     
                     DO i = 1, m
                        IF (1.0_dp-y(i) > 0.0_dp) THEN
                            grad(i) = -1.0_dp
                        END IF
                     END DO                 
                   !-------------------------------------   
                   
                   !-------------------------------------
                   !             Function S5
                   !-------------------------------------
                   CASE(5)
                     
                     DO i = 1, m
                        IF (EXP(y(i)) >= -y(i)) THEN
                            grad(i) = EXP(y(i))
                        ELSE
                            grad(i) = -1.0_dp
                        END IF
                     END DO 
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function S6
                   !-------------------------------------
                   CASE(6)
                     
                     DO i = 1, m
                        IF (EXP(-y(i)) > y(i)) THEN
                            grad(i) = -EXP(-y(i))
                        ELSE
                            grad(i) = 1.0_dp
                        END IF
                     END DO 
                   !-------------------------------------    
                   
                   !-------------------------------------
                   !             Function S7
                   !-------------------------------------
                   CASE(7)
                     
                     DO i = 1, m
                        IF (y(i) > 0.0_dp) THEN
                            grad(i) = 1.0_dp
                        ELSE    
                            grad(i) = -1.0_dp
                        END IF
                        IF (20.0_dp*(y(i)**2-y(i)-1.0_dp) > 0.0_dp) THEN
                            grad(i) = grad(i) + 40.0_dp*y(i)-20.0_dp
                        END IF
                     END DO 
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function S8
                   !-------------------------------------
                   CASE(8)
                     
                     DO i = 1, m
                        apu1 = y(i)+1.0_dp
                        apu2 = -2.0_dp*y(i)-2.0_dp
                        apu3 = 0.5_dp*y(i)+5.0_dp
                        apu = MAX(apu1, apu2, apu3)
                        IF (apu == apu1) THEN
                            grad(i) = 1.0_dp
                        ELSE IF (apu == apu2) THEN
                            grad(i) = -2.0_dp
                        ELSE
                            grad(i) = 0.5_dp
                        END IF
                     END DO
                   !-------------------------------------

                   !-------------------------------------
                   ! ************************************ 
                   !           Used in GROUP 1
                   !
                   !   i.e. block-separable problems
                   ! ************************************ 
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function B1
                   !-------------------------------------
                   CASE(1001)             
                     grad = 0.0_dp
                     ind = 1
                     
                     apu = -100.0_dp
                     DO i = 1, m
                        IF (apu < y(i)**2 ) THEN 
                            apu = y(i)**2
                            ind = i
                        END IF
                     END DO
                     grad(ind) = grad(ind) + 2.0_dp*(m+1)*y(ind)
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function B2
                   !-------------------------------------
                   CASE(1002) 
                     grad = 0.0_dp
                     ind = 1
                     
                     apu = -100.0_dp
                     DO i = 1, m
                        IF (apu < y(i)**2) THEN 
                            apu = y(i)**2
                            ind = i
                        END IF
                     END DO
                     grad(ind) = grad(ind) + 2.0_dp*y(ind)
                   !-------------------------------------              

                   !-------------------------------------
                   !             Function B3
                   !-------------------------------------
                   CASE(1003)    
                     grad = 2.0_dp * y
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function B4
                   !-------------------------------------
                   CASE(1004) 
                     grad = 0.0_dp
                     ind = 1
                     
                     apu = -100.0_dp
                     DO i = 1, m
                        IF (apu < ABS(y(i))) THEN 
                            apu = ABS(y(i))
                            ind = i
                        END IF
                     END DO
                     IF (y(ind) <= 0.0_dp) THEN 
                        grad(ind) = grad(ind) - m
                     ELSE 
                        grad(ind) = grad(ind) + m
                     END IF 
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function B5
                   !-------------------------------------
                   CASE(1005)             
                     grad = 0.0_dp
                     
                     DO i = 1, m
                        IF (y(i) <= 0.0_dp) THEN 
                            grad(i) = grad(i) - 1.0_dp
                        ELSE
                            grad(i) = grad(i) + 1.0_dp
                        END IF
                     END DO
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function B6
                   !-------------------------------------
                   CASE(1006)             
                     grad = 0.0_dp
                     
                     DO i = 2, m
                        IF ( y(i) - y(i-1) <= 0.0_dp) THEN
                            grad(i-1) = grad(i-1) + 1.0_dp 
                            grad(i) = grad(i) - 1.0_dp 
                        ELSE
                            grad(i-1) = grad(i-1) - 1.0_dp  
                            grad(i) = grad(i) + 1.0_dp     
                        END IF
                     END DO
                   !-------------------------------------   
                   
                   !-------------------------------------
                   !             Function B7
                   !-------------------------------------
                   CASE(1007)
                     grad = 0.0_dp
                     
                     apu1 = -100.0_dp
                     DO l = 1, m
                        apu2 = 0.0_dp
                        DO i = 1, m
                            apu2 = apu2 + (y(i) / REAL(l+i-1, KIND=dp))
                        END DO
                        IF (apu1 < ABS(apu2)) THEN
                            apu1 = ABS(apu2)
                            ind = l
                            f_sign1 = SIGN(1.0_dp,apu2)
                        END IF
                     END DO
                     DO i = 1, m
                        grad(i) = grad(i)+f_sign1/REAL(ind+i-1, KIND=dp)
                     END DO 
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function B8
                   !-------------------------------------
                   CASE(1008)             
                     grad = 0.0_dp
                     
                     DO i = 1, m-1
                        apu1 = -y(i)-y(i+1)
                        apu2 = -y(i)-y(i+1)+y(i)*y(i)+y(i+1)*y(i+1)-1.0_dp
                        IF (apu1 >= apu2) THEN
                            grad(i) = grad(i) - 1.0_dp
                            grad(i+1) = grad(i+1) - 1.0_dp
                        ELSE
                            grad(i) = grad(i)-1.0_dp+2.0_dp*y(i)
                            grad(i+1) = grad(i+1)-1.0_dp+2.0_dp*y(i+1)
                        END IF
                     END DO
                   !-------------------------------------

                   !-------------------------------------
                   !             Function B9
                   !-------------------------------------
                   CASE(1009)             
                     grad = 0.0_dp
                     
                     DO i = 1, m-1
                        apu1 = y(i)*y(i)*y(i)*y(i)+y(i+1)*y(i+1)
                        apu2 = (2.0_dp-y(i))*(2.0_dp-y(i))+(2.0_dp-y(i+1))*(2.0_dp-y(i+1))
                        apu3 = 2.0_dp*EXP(-y(i)+y(i+1))
                        apu = MAX(apu1, apu2, apu3)
                        IF (apu == apu1) THEN
                            grad(i) = grad(i)+4.0_dp*y(i)*y(i)*y(i)
                            grad(i+1) = grad(i+1)+2.0_dp*y(i+1)
                        ELSE IF (apu == apu2) THEN
                            grad(i) = grad(i)+2.0_dp*y(i)-4.0_dp
                            grad(i+1) = grad(i+1)+2.0_dp*y(i+1)-4.0_dp
                        ELSE
                            grad(i) = grad(i) - 2.0_dp*EXP(-y(i)+y(i+1))
                            grad(i+1) = grad(i+1) + 2.0_dp*EXP(-y(i)+y(i+1))
                        END IF   
                     END DO
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function B10
                   !-------------------------------------
                   CASE(1010)             
                     grad = 0.0_dp
                     
                     apu1 = 0.0_dp
                     apu2 = 0.0_dp
                     apu3 = 0.0_dp
                     DO i = 1, m-1
                        apu1 = apu1 + y(i)*y(i)*y(i)*y(i)+y(i+1)*y(i+1)
                        apu2 = apu2 + (2.0_dp-y(i))*(2.0_dp-y(i))+(2.0_dp-y(i+1))*(2.0_dp-y(i+1))
                        apu3 = apu3 + 2.0_dp*EXP(-y(i)+y(i+1))
                     END DO
                     apu = MAX(apu1, apu2, apu3)
                     IF (apu == apu1) THEN
                        DO i = 1, m-1
                            grad(i) = grad(i)+4.0_dp*y(i)*y(i)*y(i)
                            grad(i+1) = grad(i+1)+2.0_dp*y(i+1)
                        END DO 
                     ELSE IF (apu == apu2) THEN
                        DO i = 1, m-1
                            grad(i) = grad(i)+2.0_dp*y(i)-4.0_dp
                            grad(i+1) = grad(i+1)+2.0_dp*y(i+1)-4.0_dp
                        END DO
                     ELSE
                        DO i = 1, m-1
                            grad(i) = grad(i) - 2.0_dp*EXP(-y(i)+y(i+1))
                            grad(i+1) = grad(i+1) + 2.0_dp*EXP(-y(i)+y(i+1))
                        END DO
                     END IF
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function B11
                   !-------------------------------------
                   CASE(1011)             
                     grad = 0.0_dp
                     
                     DO i = 1, m
                        IF ( y(i) < 0.0_dp ) THEN 
                           grad(i) = grad(i)-1.0_dp
                        ELSE
                           grad(i) = grad(i)+1.0_dp
                        END IF
                        IF (20.0_dp*(y(i)*y(i)-y(i)-1.0_dp) > 0.0_dp) THEN
                           grad(i) = grad(i) + 20.0_dp*(2.0_dp*y(i)-1.0_dp)
                        END IF
                     END DO
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function B12
                   !-------------------------------------
                   CASE(1012)
                     grad = 0.0_dp
                     
                     apu1 = -100.0_dp
                     DO l = 1, m
                        apu2 = 0.0_dp
                        DO i = 1, m
                            apu2 = apu2 + (y(i) / REAL(l+i-1, KIND=dp))
                        END DO
                        IF (apu1 < ABS(apu2)) THEN
                            apu1 = ABS(apu2)
                            ind = l
                            f_sign1 = SIGN(1.0_dp,apu2)
                        END IF
                     END DO
                     DO i = 1, m
                        grad(i) = grad(i)+m*f_sign1/REAL(ind+i-1, KIND=dp)
                     END DO 
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function B13
                   !-------------------------------------
                   CASE(1013)             
                     grad = 0.0_dp
                     
                     DO l = 1, m
                        apu = 0.0_dp
                        DO i = 1, m
                            apu = apu + y(i) / REAL(l+i-1, KIND=dp)
                        END DO
                        f_sign1 = SIGN(1.0_dp,apu)
                        DO i = 1, m
                            grad(i) = grad(i) + f_sign1/REAL(l+i-1, KIND=dp)
                        END DO
                     END DO
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function B14
                   !-------------------------------------
                   CASE(1014)             
                     grad = 0.0_dp
                     
                     apu0 = -100.0_dp
                     DO i = 1, m-1
                        apu1 = y(i)*y(i)*y(i)*y(i)+y(i+1)*y(i+1)
                        apu2 = (2.0_dp-y(i))*(2.0_dp-y(i))+(2.0_dp-y(i+1))*(2.0_dp-y(i+1))
                        apu3 = 2.0_dp*EXP(-y(i)+y(i+1))
                        apu = MAX(apu1, apu2, apu3)
                        IF (apu0 < apu) THEN
                            apu0 = apu
                            ind = i
                            IF (apu == apu1) THEN
                                ind2 = 1
                            ELSE IF (apu == apu2) THEN
                                ind2 = 2
                            ELSE
                                ind2 = 3
                            END IF
                        END IF
                     END DO
                     IF (ind2 == 1) THEN
                        grad(ind) = grad(ind)+(m-1)*(4.0_dp*y(ind)*y(ind)*y(ind))
                        grad(ind+1) = grad(ind+1)+(m-1)*2.0_dp*y(ind+1)
                     ELSE IF (ind2 == 2) THEN
                        grad(ind) = grad(ind)+(m-1)*(2.0_dp*y(ind)-4.0_dp)
                        grad(ind+1) = grad(ind+1)+(m-1)*(2.0_dp*y(ind+1)-4.0_dp)
                     ELSE
                        grad(ind) = grad(ind) - (m-1)*(2.0_dp*EXP(-y(ind)+y(ind+1)))
                        grad(ind+1) = grad(ind+1) + (m-1)*(2.0_dp*EXP(-y(ind)+y(ind+1)))
                     END IF
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function B15
                   !-------------------------------------
                   CASE(1015)             
                     grad = 0.0_dp
                     
                     apu = 0.0_dp
                     DO i = 1, m-1
                        apu = apu + y(i)**2 + (y(i+1)-1.0_dp)**2
                        apu = apu + y(i+1) - 1.0_dp
                     END DO
                     apu = 2.0_dp * apu
                     IF (apu > 0.0_dp) THEN
                        DO i = 1, m-1
                            grad(i) = grad(i) + 4.0_dp *y(i)
                            grad(i+1) = grad(i+1) + 4.0_dp *(y(i+1)-1.0_dp) +2.0_dp
                        END DO
                     END IF 
                   !-------------------------------------
                   
                   !-------------------------------------
                   !             Function B16
                   !-------------------------------------
                   CASE(1016)             
                     grad = 0.0_dp
                     
                     DO i = 1, m-1
                        grad(i) = grad(i) + 2.0_dp * y(i) 
                        grad(i+1) = grad(i+1) + 2.0_dp * (y(i+1) -1.0_dp) +1.0_dp
                     END DO
                   !-------------------------------------    
                        
                END SELECT  

           END FUNCTION subgradient_f2_block
           
           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
           
           
           FUNCTION summa(y, j, user_n) RESULT(f)           
                !
                ! Calculates the sum used in f_1 and f_2 for parameter t_j
                !
                ! NOTICE: The dimension of 'y' has to be 'n'. 'j' needs to be an integer from interval [1,20]
                !
                IMPLICIT NONE
                !**************************** NEEDED FROM USER *************************************
                REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: y    ! a point where the function value of the DC component f_2 is calculated
                INTEGER, INTENT(IN) ::  j    ! used to determine the parameter t_j
                INTEGER, INTENT(IN) :: user_n                   ! the dimension of the problem              
                !**************************** OTHER VARIABLES **************************************
                REAL(KIND=dp) :: f, t, apu                          ! the function value of the sum used in f_1 and parameter t_j
                INTEGER :: i                                    ! help variable
                
                f = 0.0_dp
                apu = 1.0_dp/ user_n
                
                DO i = 1, user_n 
                    t = (0.05_dp * j )**(i-1)
                    f = f + (y(i)-apu)*t            
                END DO 
            
           END FUNCTION summa

    
      END MODULE functions     