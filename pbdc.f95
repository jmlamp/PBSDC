      MODULE pbdc
      
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !| |                                                                                  | |
        !| |             THE PROXIMAL BUNDLE METHOD FOR NONSMOOTH DC OPTIMIZATION             | | 
        !| |                                 (version 2.5)                                    | |
        !| |                                                                                  | |
        !| |             NEW features :                                                       | |
        !| |                                                                                  | |
        !| |                * Possibility to use simple stepsize determination after          | |
        !| |                  each 'main iteration'                                           | |
        !| |                                                                                  | |
        !| |                * During each round of 'main iteration' it is possible to         | |
        !| |                  utilize OpenMP to calculate subproblems in parallel.            | |
        !| |                  In order to do this, you need to use '-fopenmp' in Makefile.    | |
        !| |                                                                                  | |
        !| |                                                                                  | |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*     
        !|                                                                                      |
        !|    The code utilizes new version of PLQDF1 by Ladislav Luksan as a quadratic solver. |
        !|                                                                                      |
        !|    The code utilizes PVMM by Ladislav Luksan as a norm minimization solver. This     | 
        !|    subroutine uses PQSUBS and MQSUBS by Ladislav Luksan. PVMM are a VARIABLE METRIC  |
        !|    ALGORITHM for UNCONSTRAINED and LINEARLY CONSTRAINED OPTIMIZATION.                |
        !|                                                                                      |
        !|    The subroutine PVMM together with PQSUBS and MQSUBS is licensed by                |
        !|    the GNU Lesser General Public License (LGPL).                                     |
        !|                                                                                      |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        
        !$USE omp_lib               ! Contains the Open MP Fortran runtime library routines (Needed in parallellization)
        
        USE constants, ONLY : dp    ! Double precision (i.e. accuracy)
        USE bundle1                 ! The BUNDLE of the DC component f_1
        USE bundle2                 ! The BUNDLE of the DC component f_2
        USE functions               ! Contains INFORMATION from the USER
        USE norm_min                ! The SOLVER for the norm minimization problem
        
        IMPLICIT NONE    
        
        EXTERNAL PLQDF1             ! The QUADRATIC SOLVER by Ladislav Luksan
        
        
        CONTAINS
                
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !| |                                                                                  | |
        !| |                            CONTAINS SUBROUTINES:                                 | | 
        !| |                                                                                  | |
        !| |        PBDC_algorithm       : The proximal bundle method for DC optimization.    | | 
        !| |        PBDC_main_iteration  : The main iteration algorithm needed in the         | |
        !| |                                proximal bundle method for DC optimization.       | |
        !| |        subproblem_solver    : The solver for subproblems in the search direction | |
        !| |                                problem. Needed in the main iteration algorithm.  | |
        !| |                                                                                  | |       
        !| |                                                                                  | |       
        !| |                            CONTAINS FUNCTIONS:                                   | |
        !| |                                                                                  | |
        !| |        select_value_t : Selects the parameter t from interval [t_min,t_max].     | |               
        !| |                                                                                  | |       
        !| |                                                                                  | |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*   
        
    
    
        
        
        ! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  
        !| | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |
        ! START ** START ** START ** START ** START ** START ** START ** START ** START ** START  
        !|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|        
        ! <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>
        !***************************************************************************************
        !  ----------------------------------------------------------------------------------  |
        !  |                                                                                |  |
        !  |          THE PROXIMAL BUNBLE ALGORITHM FOR NONSMOOTH DC OPTIMIZATION           |  |
        !  |                                                                                |  |
        !  ----------------------------------------------------------------------------------  |
        !***************************************************************************************
        
           SUBROUTINE PBDC_algorithm( f_solution, mit, CPUtime, time_limit, problem1, problem2, &
                            &  user_n, mrounds, termination, counter, agg_used, stepsize_used, iprint, &
                            & starting_points, start_ind, nochange_mainit, nochange_tol )
            ! 
            ! Solves the unconstrained nonsmooth DC minimization problem
            !
            ! INPUT:  * 'user_n'          : the dimension of the problem
            !         * 'problem1'        : the f1 of the problem   
            !         * 'problem2'        : the f2 of the problem   
            !         * 'mit'             : The maximum number of 'main iterations'        
            !         * 'mrouds'          : The maximum number of rounds during one 'main iteration'
            !         * 'agg_used'        : If .TRUE. then aggregation is used in the algorithm            
            !         * 'stepsize_used'   : If .TRUE. then simple stepsize determination is used in the algorithm after each 'main iteration'          
            !         * 'iprint'          : Specifies the print 
            !         * 'starting_points' : Table containing different starting points          
            !         * 'start_ind'       : The index of the used starting point         
            !         * 'time_limit'      : CPU time limit for execution    
            !         * 'nochange_mainit' : 0=not used, integer value 'k' = algorithm is stopped if after 'k' main iterations change in the objective value is less than 'nochange_tol'      
            !         * 'nochange_tol'    : tolerance for change in the objective value              
            !
            ! OUTPUT: * 'x_solution'  : The solution obtained to the minimizatin problem
            !         * 'f_solution'  : The objective function value at the solution 'x_solution'
            !         * 'termination' : The cause of the termination in the algorithm
            !         * 'counter'     : Gives the values of different counters
            !         * 'CPUtime'     : CPU time
            !
            ! NOTICE: * The dimension of the vector 'counter' has to be 6.
            !         * 'mit' and 'mrounds' has to be integers.
            !         * IF ('mit' <= 0) THEN DEFAULT value 1000 is used
            !         * IF ('mrounds' <= 0 ) THEN DEFAULT value 5000 is used.
            !         * 'iprint' has to be -4, -3, -2, -1, 0, 1, 2, 3 or 4 (5). If it is NOT then DEFAULT value 1 is used. 
            !         * If 'iprint = 5' then everything is printed step by step (this is NOT recommended)
            !
            !***********************************************************************************
               IMPLICIT NONE
            !**************************** NEEDED FROM USER ************************************* 
               
               REAL(KIND=dp), DIMENSION(:,:), INTENT(IN) :: starting_points  ! the table containing different starting points
               INTEGER, INTENT(IN)                       :: start_ind        ! the index of the used starting point
               
               REAL(KIND=dp), INTENT(OUT) :: f_solution                ! the objective function value at the solution 'x_solution'
               REAL(KIND=dp), INTENT(OUT) :: CPUtime                   ! the CPU time
               REAL(KIND=dp), INTENT(IN) :: time_limit                 ! the CPU time limit for execution
               
               INTEGER, INTENT(IN) :: user_n                     ! the dimension
               INTEGER, DIMENSION(:), INTENT(IN) :: problem1     ! the f1
               INTEGER, DIMENSION(:), INTENT(IN) :: problem2     ! the f2
               
               INTEGER, INTENT(INOUT) :: mit                  ! the maximum number of 'main iterations'
               INTEGER, INTENT(INOUT) :: mrounds              ! the maximum number of rounds during one 'main iteration'
               
               INTEGER, INTENT(IN) :: nochange_mainit       ! 0=not used, integer value 'k' = algorithm is stopped if after 'k' main iterations change in the objective value is less than 'nochange_tol'
               REAL(KIND=dp), INTENT(IN) :: nochange_tol    ! tolerance for change in the objective value
               
               INTEGER, INTENT(OUT) :: termination        ! 1 - the stopping condition satisfied (i.e. criticality)
                                                          ! 2 - the approximate stopping condition satisfied (i.e. eps-criticality)
                                                          ! 3 - the maximum number 'mrounds' of rounds executed in one main iteration
                                                          ! 4 - the maximum number of 'main iterations' is executed  
                                                          ! 5 - the executed CPU time exceeded the CPU time limit       
                                                          ! 7 - after 'nochange_mainit' main iterations change in the objective value is less than 'nochange_tol'                                                         
                                                          
               INTEGER, DIMENSION(8), INTENT(OUT) :: counter  ! contains the values of different counteres: 
                                                              !   counter(1) = iter_counter         the number of 'main iterations' executed
                                                              !   counter(2) = subprob_counter      the number of subproblems solved
                                                              !   counter(3) = f_counter            the number of function values evaluated for DC component
                                                              !   counter(4) = subgrad1_counter     the number of subgradients calculated for f_1 
                                                              !   counter(5) = subgrad2_counter     the number of subgradients calculated for f_2 
                                                              !   counter(6) = stop_cond_counter    the number of times approximate stopping condition was tested during the algorithm 

               LOGICAL, INTENT(IN) :: agg_used          ! .TRUE. if aggregation is used in the algorithm. Otherwise .FALSE.                                                           
               LOGICAL, INTENT(IN) :: stepsize_used     ! .TRUE. if simple stepsize determination is used in the algorithm. Otherwise .FALSE.                                                             
            
               INTEGER, INTENT(INOUT) :: iprint ! variable that specifies print option:
                                                !   iprint = 0 : print is suppressed
                                                !   iprint = 1 : basic print of final result 
                                                !   iprint = -1: basic print of final result (without the solution vector)
                                                !   iprint = 2 : extended print of final result 
                                                !   iprint = -2: extended print of final result (without the solution vector)
                                                !   iprint = 3 : basic print of intermediate results and extended print of final results
                                                !   iprint = -3: basic print of intermediate results and extended print of final results (without the solution vector)
                                                !   iprint = 4 : extended print of intermediate results and extended print of final results 
                                                !   iprint = -4: extended print of intermediate results and extended print of final results (without the solution vectors)
                                                !   iprint = 5 : prints each step of the bundle algorithm (i.e. everything is printed step by step) 
           
           
           !***************************** LOCAL VARIABLES ************************************  
           
               TYPE(kimppu1) :: B1             ! The bundle B_1 for the DC component f_1
               TYPE(kimppu2) :: B2             ! The bundle B_2 for the DC component f_2

               REAL(KIND=dp), DIMENSION(user_n) :: x_current     ! the current iteration point (the dimension 'user_n' is the number of variables)
               REAL(KIND=dp), DIMENSION(user_n) :: x_new         ! the new iteration point (obtained from the previous 'main iteration')  
               
               REAL(KIND=dp), DIMENSION(user_n) :: x_0           ! the starting point
               REAL(KIND=dp), DIMENSION(user_n) :: x_solution    ! the solution obtained to the problem
           
               REAL(KIND=dp), DIMENSION(user_n) :: grad1, grad2  ! subgradients (the dimenison 'user_n' is the length of subgradient)       
               REAL(KIND=dp), DIMENSION(user_n) :: dd            ! the search direction obtained from the 'main iteration': dd=x_new -x_current
                                                                 ! (the dimension 'user_n' is the length of subgradient)               

               REAL(KIND=dp) :: L1, L2         ! Lipschitz constants of DC components f_1 and f_2
               REAL(KIND=dp) :: crit_tol, eps  ! 'crit_tol'=the criticality tolerance and 'eps'=the proximity measure (needed in stopping conditions)
               REAL(KIND=dp) :: m              ! the descent parameter
               REAL(KIND=dp) :: r_dec, r_inc   ! 'r_dec'=the decrease parameter and 'R_inc'=the increase parameter         

               REAL(KIND=dp) :: c              ! 'c'=the decrease parameter
               
               REAL(KIND=dp) :: f_0                         ! the value of the objective function f=f_1-f_2 at the starting point x_0
               REAL(KIND=dp) :: f1_current, f2_current      ! the value of f_1 and f_2 at the current solution  
               REAL(KIND=dp) :: f1_new, f2_new              ! the value of f_1 and f_2 at the new iteration point x_new (obtained from the previous main iteration)        
               REAL(KIND=dp) :: change                      ! 'change' = f(x_new) - f(x_current)  (i.e. the change in the objective function value)
               REAL(KIND=dp) :: change1                     ! 'change1' = f1(x_new) - f1(x_current)  (i.e. the change in the DC component f1)
               REAL(KIND=dp) :: change2                     ! 'change2' = f2(x_new) - f2(x_current)  (i.e. the change in the DC component f2)
    
               REAL(KIND=dp) :: start_time, finish_time                   ! start and finish CPU time
               REAL(KIND=dp) :: start_time_main_it, finish_time_main_it   ! start and finish CPU time in one 'main iteration'
               
               REAL(KIND=dp) :: middle_time                               ! Middle CPU time     

               REAL(KIND=dp), DIMENSION(nochange_mainit) :: fk_table  ! table for saving objective function values at the previous iteration points            
               
               REAL(KIND=dp) :: elapsed_time                  ! elapsed 'clock' time in seconds
               INTEGER :: clock_start, clock_end, clock_rate  ! start and finish 'clock' time   

               INTEGER :: size_b1 , size_b2    ! The biggest possible size of the bundles B_1 and B_2 
               
               INTEGER  :: iter_counter        ! the number of main iterations executed 
               INTEGER  :: subprob_counter     ! the number of subproblems solved during the execution of the bundle algorithm
               INTEGER  :: stop_cond_counter   ! the number of times approximate stopping condition is tested during the algorithm
               INTEGER  :: f_counter           ! the number of function values evaluated for a DC component during 
                                               ! the execution of the bundle algorithm (same for the DC component f_1 and f_2)
                                               
               INTEGER  :: subgrad1_counter    ! the number of subgradients calculated for f_1 during the bundle algorithm
               INTEGER  :: subgrad2_counter    ! the number of subgradients calculated for f_2 during the bundle algorithm                 
               
               INTEGER :: help_mit_iter_counter    ! the number of iteration rounds executed in one 'main iteration'
               INTEGER :: help_subprob_counter     ! the number of subproblems solved in one 'main iteration'
               INTEGER :: help_f_counter           ! the number of function values evaluated for a DC component in one 'main iteration' (same for f_1 and f_2)
               INTEGER :: help_subgrad1_counter    ! the number of subgradients calculated for f_1 in one 'main iteration'
               INTEGER :: help_subgrad2_counter    ! the number of subgradients calculated for f_2 in one 'main iteration'
               INTEGER :: help_stop_cond_counter   ! the number of times approximate stopping condition was tested during the 'main iteration'
               
               INTEGER :: reason_for_stop       ! the reason for stop during the 'main iteration'
                                                ! 0 - a new iteration point found
                                                ! 1 - the stopping condition satisfied (criticality)
                                                ! 2 - the approximate stopping condition satisfied (eps-criticality)
                                                ! 3 - the biggest possible number of rounds executed in the 'main iteration'
               
               INTEGER :: max_threads           ! the maximum number of threads that can be used in parallellization
               INTEGER :: threads               ! the number of threads used in parallellization
               INTEGER :: max_sub_prob          ! the maximum number of subproblems solved at the same time
               
               INTEGER :: i, j                  ! help variables    
                       
               LOGICAL :: stop_alg      ! .TRUE. if the proximal bundle algorithm can be terminated                                                        
               
               CALL cpu_time(start_time)                ! Start CPU timing     
               CALL cpu_time(middle_time)                ! Start CPU timing     

               CALL SYSTEM_CLOCK(COUNT_RATE=clock_rate) ! Find the rate
               CALL SYSTEM_CLOCK(COUNT=clock_start)     ! Start timing 'clock' time            


!--------------------------------------------------------------------------------------
          
! The starting point for the test problem is selected from the 'starting_points' table
                
                ! The starting point of the first block is multiplied for other blocks
                DO i = 1,number_of_blocks
                   DO j = 1,size_of_block(1)
                       x_0(size_of_block(1)*(i-1)+j) = starting_points(start_ind,j)
                   END DO
                END DO
                
!--------------------------------------------------------------------------------------
               
           !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^     
           ! <>  <>  <>  <>  <>  <>  <>  MAIN ITERATION STARTS  <>  <>  <>  <>  <>  <>  <>  <>
           !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^      
               
               !************ PARAMETER VALUES FROM USER ARE LOOKED THROUGH ********************
               !
               !     IF VALUES ARE SET WRONG THEN DEFAULT PARAMETER VALUES ARE USED
               !
               ! ***** Lipschitz constant 'L1' *****
               IF (user_L1 <= 0.0_dp) THEN
                   L1 = 1000.0_dp   
               ELSE    
                   L1 = user_L1
               END IF
              
               ! ***** Lipschitz constant 'L2' *****
               IF (user_L2 <= 0.0_dp) THEN
                   L2 = 1000.0_dp
               ELSE
                   L2 = user_L2
               END IF
              
               ! ***** descent parameter 'm' *****
               IF ( (user_m<=0.0_dp) .OR. (user_m>=1.0_dp) ) THEN
                   m = 0.2_dp
               ELSE       
                   m = user_m 
               END IF
              
               ! ***** criticality tolerance 'crit_tol' *****
               IF (user_crit_tol <= 0.0_dp) THEN
               
                   IF (user_n < 150) THEN 
                      crit_tol = user_n*0.005_dp
                   ELSE IF ( (150 <= user_n) .AND. (user_n <=200) ) THEN 
                      crit_tol = user_n*0.015_dp
                   ELSE IF (user_n > 200) THEN 
                      crit_tol = user_n*0.05_dp
                   END IF
                   
               ELSE 
                   crit_tol = user_crit_tol
               END IF
              
               ! ***** proximity parameter 'eps' *****
               IF (user_eps <= 0.0_dp) THEN
                   eps = 0.1_dp
               ELSE
                   eps = user_eps
               END IF
              
               ! ***** decrease parameter 'r_dec' *****
               IF ((user_r_dec <= 0.0_dp) .OR. (user_r_dec >= 1.0_dp) ) THEN
                   
                   IF( user_n < 10) THEN 
                      r_dec = 0.75_dp
                   ELSE IF ( user_n == 10) THEN 
                      r_dec = 0.66_dp
                   ELSE IF ( user_n == 11) THEN 
                      r_dec = 0.68_dp   
                   ELSE IF ( user_n == 12) THEN 
                      r_dec = 0.70_dp
                   ELSE IF ( user_n == 13) THEN 
                      r_dec = 0.72_dp
                   ELSE IF ( user_n == 14) THEN 
                      r_dec = 0.73_dp
                   ELSE IF ( user_n == 15) THEN 
                      r_dec = 0.75_dp
                   ELSE IF ( user_n == 16) THEN 
                      r_dec = 0.76_dp
                   ELSE IF ( user_n == 17) THEN 
                      r_dec = 0.77_dp
                   ELSE IF ( user_n == 18) THEN 
                      r_dec = 0.78_dp
                   ELSE IF ( user_n == 19) THEN 
                      r_dec = 0.79_dp                     
                   ELSE IF ( user_n == 20 .OR. user_n == 21 ) THEN 
                      r_dec = 0.80_dp    
                   ELSE IF ( user_n == 22 ) THEN 
                      r_dec = 0.81_dp    
                   ELSE IF ( user_n == 23 .OR. user_n == 24 ) THEN 
                      r_dec = 0.82_dp    
                   ELSE IF ( user_n == 25 .OR. user_n == 26 ) THEN 
                      r_dec = 0.83_dp           
                   ELSE IF ( user_n == 27 .OR. user_n == 28 ) THEN 
                      r_dec = 0.84_dp                         
                   ELSE IF ( user_n == 29 .OR. user_n == 30 ) THEN 
                      r_dec = 0.85_dp                   
                   ELSE IF ( user_n >= 31 .AND. user_n <= 33 ) THEN 
                      r_dec = 0.86_dp       
                   ELSE IF ( user_n >= 34 .AND. user_n <= 36 ) THEN 
                      r_dec = 0.87_dp  
                   ELSE IF ( user_n >= 37 .AND. user_n <= 40 ) THEN 
                      r_dec = 0.88_dp       
                   ELSE IF ( user_n >= 41 .AND. user_n <= 44 ) THEN 
                      r_dec = 0.89_dp               
                   ELSE IF ( user_n >= 45 .AND. user_n <= 50 ) THEN 
                      r_dec = 0.90_dp                             
                   ELSE IF ( user_n >= 51 .AND. user_n <= 57 ) THEN 
                      r_dec = 0.91_dp  
                   ELSE IF ( user_n >= 58 .AND. user_n <= 66 ) THEN 
                      r_dec = 0.92_dp   
                   ELSE IF ( user_n >= 67 .AND. user_n <= 78 ) THEN 
                      r_dec = 0.93_dp                 
                   ELSE IF ( user_n >= 79 .AND. user_n <= 94 ) THEN 
                      r_dec = 0.94_dp                     
                   ELSE IF ( user_n >= 95 .AND. user_n <= 119 ) THEN 
                      r_dec = 0.95_dp
                   ELSE IF ( user_n >= 120 .AND. user_n <= 161 ) THEN 
                      r_dec = 0.96_dp
                   ELSE IF ( user_n >= 162 .AND. user_n <= 244 ) THEN 
                      r_dec = 0.97
                   ELSE IF ( user_n >= 245 .AND. user_n <= 299 ) THEN 
                      r_dec = 0.98
                   ELSE IF ( user_n >= 300 ) THEN 
                      r_dec = 0.99
                   END IF 
                   
               ELSE
                   r_dec = user_r_dec
               END IF

               ! ***** increase parameter 'r_inc' *****
               IF ( user_r_inc <= 1.0_dp ) THEN
                   r_inc = 10000000.0_dp 
               ELSE
                   r_inc = user_r_inc
               END IF

               ! ***** bundle B1 size *****
               IF ( user_size_b1 <= 0 ) THEN
                   size_b1 = MIN(user_n+5,100) 
               ELSE
                   size_b1 = user_size_b1
               END IF

               ! ***** bundle B2 size *****
               IF ( user_size_b2 <= 0 ) THEN
                   size_b2 = 3
               ELSE
                   size_b2 = user_size_b2
               END IF          

               ! ***** decrease parameter 'c' *****               
               IF ((user_c<0.0_dp) .OR. (user_c>1.0_dp)) THEN
                   c = 0.1_dp
               ELSE
                   c = user_c
               END IF              
               
               !----------------------------------------------------------------------------
               ! ***** maximum number of main iterations 'mit' *****
               IF ( mit <= 0 ) THEN
                   mit = 1000
               END IF   
               
               ! ***** maximum number of rounds 'mrounds' in the 'main iteration' *****
               IF ( mrounds <= 0 ) THEN
                   mrounds = 5000
               END IF   
               
               ! ***** print option 'iprint' *****
               IF ( (iprint < -4 ) .OR. (iprint > 5) ) THEN   !executed if print value is wrong
                   iprint = 1
               END IF              
              !----------------------------------------------------------------------------
              
           !_______________________________________________________________________________
           !************************ STEP 0: PARAMETER INITIALIZATION *********************                
           
               x_current = x_0         ! the current iteration point is the starting point x_0
               f1_current = f1(x_0,problem1,user_n)    ! the value of the DC component f_1 at x_0
               f2_current = f2(x_0,problem2,user_n)    ! the value of the DC component f_2 at x_0
               f_counter = 1           ! one function value evaluated for both f_1 and f_2
           
               f_0 = f1_current - f2_current  ! the value of the objective function at the starting point x_0 
               
               grad1 = subgradient_f1(x_0,problem1,user_n)    ! the subgradient of f_1 at x_0
               subgrad1_counter = 1           ! one subgradient calculated for f_1
               grad2 = subgradient_f2(x_0,problem2,user_n)    ! the subgradient of f_2 at x_0
               subgrad2_counter = 1           ! one subgradient calculated for f_1
               
               ! The bundles B_1 and B_2 are initialized
               CALL init_bundle_b1(B1, size_b1, user_n)    
               CALL init_bundle_b2(B2, size_b2, user_n)
               
               ! The first bundle element is added into the bundle B_1 and B_2 (i.e. the one corresponding to the starting point)
               CALL add_first_element_b1(B1, grad1)
               CALL add_first_element_b2(B2, grad2)
               
               iter_counter = 0             ! the number of 'main iterations' executed so far is zero
               subprob_counter = 0          ! the number of 'subproblems' solved so far is also zero
               stop_cond_counter = 0        ! the number of times approximate stopping condition is tested is zero
               stop_alg = .FALSE.           ! we cannot stop the proximal bundle algorithm
     
! --- --- --- Needed in OpenMP when we use PARALLELLIZATION --- --- ---      
               !$max_threads = omp_get_max_threads()
               !$max_sub_prob = give_max_size_b2(B2)+1
               !$threads = MIN(max_threads, max_sub_prob)
               !$CALL omp_set_num_threads(threads)              
! ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----            
                   
  
              IF ( (ABS(iprint) >= 3) .AND. (ABS(iprint) <= 5) ) THEN        ! the basic/extended print of indermediate results
                        
                   WRITE(*,*) 'Main Iter:', iter_counter, 'f(x):', f1_current-f2_current    
              END IF       
                   
           !_______________________________________________________________________________
           !************************ STEP 0: END ******************************************    
           
           !----------------------------------------------------------------------------------   
           ! <>  <>  <>  <>  <>  <>  <> REPEATABLE PART BEGINS <>  <>  <>  <>  <>  <>  <>  <>
           !----------------------------------------------------------------------------------  
           
            DO WHILE ( ( .NOT. stop_alg) .AND. (iter_counter < mit ) .AND.  &
                          ((middle_time-start_time) < time_limit) )  ! is repeated until the bundle algorithm can be terminated 
                                                                         
               !_______________________________________________________________________________
               !************************ STEP 1: MAIN ITERATION ******************************* 

                iter_counter = iter_counter + 1                     ! a new 'main iteration' is executed
                
                IF(iprint==5) THEN  !prints everything
                    WRITE(*,*) '-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'
                    WRITE(*,*) '         MAIN ITERATION', iter_counter, ' STARTS'
                    WRITE(*,*) '-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'
                    WRITE(*,*) ' '
                END IF 

                CALL cpu_time(start_time_main_it)
                
                CALL PBDC_main_iteration( x_current, f1_current, f2_current,&
                            &  x_new, f1_new, f2_new, f_0, B1, B2,  &
                            & L1, L2, crit_tol, eps, m, c, r_dec, r_inc, & 
                            & mrounds, iprint, reason_for_stop, help_mit_iter_counter, &
                            & help_subprob_counter, help_f_counter, help_subgrad1_counter, &
                            & help_subgrad2_counter, help_stop_cond_counter, agg_used, &
                            & stepsize_used, problem1, problem2, user_n)    

                IF(iprint==5) THEN  !prints everything
                    WRITE(*,*) '-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'
                    WRITE(*,*) '         MAIN ITERATION', iter_counter, ' ENDS'
                    WRITE(*,*) '-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'
                    WRITE(*,*) ' '
                END IF 
                        
                CALL cpu_time(finish_time_main_it)
                
                ! Different counters are updated
                subprob_counter = subprob_counter + help_subprob_counter
                f_counter = f_counter + help_f_counter
                stop_cond_counter = stop_cond_counter + help_stop_cond_counter
                subgrad1_counter = subgrad1_counter + help_subgrad1_counter
                subgrad2_counter = subgrad2_counter + help_subgrad2_counter 

            
               !_______________________________________________________________________________
               !************************ STEP 1: END ****************************************** 
               
               IF (reason_for_stop == 0) THEN   ! in this case, a new iteration point is found in the 'main iteration'
               !_______________________________________________________________________________
               !************************ STEP 2: BUNDLE UPDATE*******************************
                
                    ! Subgradients of f_1 and f_2 at the new point x_new
                    grad1 = subgradient_f1(x_new,problem1,user_n)
                    subgrad1_counter = subgrad1_counter + 1
                    change1 = f1_new - f1_current
                    
                    grad2 = subgradient_f2(x_new,problem2,user_n)
                    subgrad2_counter = subgrad2_counter + 1
                    change2 = f2_new - f2_current
                    
                    dd = x_new - x_current                              ! the search direction
                    change = f1_new - f2_new - (f1_current -f2_current) ! the change in the objective function value
                
                    CALL update_b1(B1, grad1, dd, change1)              ! bundle update for B_1
                    CALL update_b2(B2, grad2, dd, change2)              ! bundle update for B_2
            
                    x_current = x_new                                   ! update of the current iteration point
                    f1_current = f1_new                                 ! update of the function value f_1
                    f2_current = f2_new                                 ! update of the function value f_2
                    
                    ! Check if the objective function value has not changed in 'nochange_mainit' iterations
                    IF (nochange_mainit > 0) THEN
                        IF (iter_counter < nochange_mainit+1) THEN
                            fk_table(iter_counter)=f1_current-f2_current
                        ELSE
                            IF (ABS(fk_table(1)-(f1_current-f2_current)) < nochange_tol) THEN
                                reason_for_stop = 7
                                stop_alg = .TRUE.   ! the minimization algorithm can be STOPPED
                            ELSE
                                DO i=1,nochange_mainit-1
                                    fk_table(i)=fk_table(i+1)
                                END DO
                                fk_table(nochange_mainit)=f1_current-f2_current
                            END IF
                        END IF
                    END IF

                    CALL cpu_time(middle_time)         ! Middle CPU timing
           
                    IF ( (ABS(iprint) >= 3) .AND. (ABS(iprint) <= 5) ) THEN        ! the basic/extended print of indermediate results
                        
                        WRITE(*,*) 'Main Iter:', iter_counter, 'f(x):', f1_current-f2_current, '  change:', change
                        IF (iprint > 3) THEN    ! IF iprint =  4 or 5 then the iteration point is printes
                            WRITE(*,*) 'The new iteration point:'
                            DO i = 1, user_n
                                WRITE(*,*) 'x(', i ,')*=',  x_current(i)
                            END DO
                            WRITE(*,*) ' '
                        END IF 
                        
                        IF ( (ABS(iprint) == 4 ) .OR. (ABS(iprint) == 5 ) ) THEN
                            WRITE(*,*) '--------------------------------------------------------------------------- '
                            WRITE(*,*) ' * * * * * DETAILED INFORMATION ABOUT THE MAIN ITERATION', iter_counter ,' * * * * * ' 
                            WRITE(*,*) 'the number of rounds needed:', help_mit_iter_counter
                            WRITE(*,*) 'the number of subproblems solved:', help_subprob_counter
                            WRITE(*,*) 'the number of function values calculated for f_1:', help_f_counter
                            WRITE(*,*) 'the number of function values calculated for f_2:', help_f_counter
                            WRITE(*,*) 'the number of subgradients calculated for f_1:', help_subgrad1_counter
                            WRITE(*,*) 'the number of subgradients calculated for f_2:', help_subgrad2_counter
                            WRITE(*,*) 'the number of times approximate stopping condition was tested:',help_stop_cond_counter 
                            WRITE(*,*) '--------------------------------------------------------------------------- '   
                            WRITE(*,*) 'CPU time used:', finish_time_main_it-start_time_main_it                             
                            WRITE(*,*) '--------------------------------------------------------------------------- '   
                            WRITE(*,*) ' '                          
                        END IF
                        

                        
                    END IF  
                            
               !_______________________________________________________________________________
               !************************ STEP 2: END ****************************************       
               
               ELSE                     ! one of the stopping conditions is fulfilled
               
                    stop_alg = .TRUE.   ! the minimization algorithm can be STOPPED
                    
                    IF((ABS(iprint) >= 4) .AND. (ABS(iprint) <= 5)) THEN ! the extended print of indermediate results
                        WRITE(*,*) ' '
                        IF ( (ABS(iprint) == 4 ) .OR. (ABS(iprint) == 5 ) ) THEN
                        WRITE(*,*) '--------------------------------------------------------------------------- '
                        WRITE(*,*) ' * * * * * DETAILED INFORMATION ABOUT THE MAIN ITERATION', iter_counter ,' * * * * * '
                        WRITE(*,*) 'the number of rounds needed:', help_mit_iter_counter
                        WRITE(*,*) 'the number of subproblems solved:', help_subprob_counter
                        WRITE(*,*) 'the number of function values calculated for f_1:', help_f_counter
                        WRITE(*,*) 'the number of function values calculated for f_2:', help_f_counter
                        WRITE(*,*) 'the number of subgradients calculated for f_1:', help_subgrad1_counter
                        WRITE(*,*) 'the number of subgradients calculated for f_2:', help_subgrad2_counter
                        WRITE(*,*) 'the number of times approximate stopping condition was tested:',help_stop_cond_counter 
                        WRITE(*,*) '--------------------------------------------------------------------------- '
                        WRITE(*,*) 'CPU time used:', finish_time_main_it-start_time_main_it                    
                        WRITE(*,*) '--------------------------------------------------------------------------- '
                        END IF      
                        WRITE(*,*) ' '
                        WRITE(*,*) 'During Main iteration round:', iter_counter
                        WRITE(*,*) 'Some of the stopping conditions is fulfilled and the algorithm is stopped.'                     
                        WRITE(*,*) ' '
                    END IF 
                    
               END IF
               
           END DO
           !----------------------------------------------------------------------------------   
           ! <>  <>  <>  <>  <>  <>  <>  REPEATABLE PART ENDS  <>  <>  <>  <>  <>  <>  <>  <> 
           !----------------------------------------------------------------------------------         
           

           termination = reason_for_stop   ! the cause of the termination
           
           IF ( iter_counter >= mit ) THEN  ! the maximum number of 'main iterations' is executed and this causes the termination
               termination = 4 
           END IF

           IF ( middle_time-start_time > time_limit ) THEN  ! the maximum number of 'main iterations' is executed and this causes the termination
               termination = 5 
           END IF           
           
           x_solution = x_current               ! the solution to the minimization problem
           f_solution = f1_current - f2_current ! the objective function value at the solution
           
           ! the values of different counters
           counter(1) = iter_counter 
           counter(2) = subprob_counter
           counter(3) = f_counter
           counter(4) = subgrad1_counter
           counter(5) = subgrad2_counter 
           counter(6) = stop_cond_counter
           counter(7) = 0
           counter(8) = 0
            
           CALL cpu_time(finish_time)         ! Stop CPU timing
           CALL SYSTEM_CLOCK(COUNT=clock_end) ! Stop timing 'clock' time
           
          ! Calculate the elapsed 'clock' time in seconds:
          elapsed_time=(1.0_dp*clock_end-clock_start)/clock_rate           

          CPUtime = finish_time-start_time 
          
           IF ( (ABS(iprint) == 1) ) THEN ! basic print of the final result
               WRITE(*,*) '---------------------------------------------------------------------------'        
               WRITE(*,*) '           * * * * * BASIC PRINT OF THE FINAL SOLUTION * * * * * '
               WRITE(*,*) '---------------------------------------------------------------------------'            
               WRITE(*,*) 'The cause of termination: ', termination
               IF (iprint > 0) THEN
                   WRITE(*,*) '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - '
                   DO i = 1, user_n
                       WRITE(*,*) 'x(', i ,')*=',  x_solution(i)
                   END DO
               END IF
               WRITE(*,*) '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - '               
               WRITE(*,*) 'f*=', f_solution     
               WRITE(*,*) '---------------------------------------------------------------------------' 
               WRITE(*,*) 'CPU time used:', finish_time-start_time                              
               WRITE(*,*) '---------------------------------------------------------------------------'
               WRITE(*,*) 'Elapsed time:', elapsed_time
               WRITE(*,*) '---------------------------------------------------------------------------'            
           END IF
           
           IF ((ABS(iprint) /= 1) .AND. (iprint /= 0)) THEN  ! extended print of the final result
               WRITE(*,*) '---------------------------------------------------------------------------'        
               WRITE(*,*) '        * * * * * EXTENDED PRINT OF THE FINAL SOLUTION * * * * * '
               WRITE(*,*) '---------------------------------------------------------------------------'            
               WRITE(*,*) 'The cause of termination: ', termination
               IF (iprint > 0) THEN
                   WRITE(*,*) '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - '
                   DO i = 1, user_n
                       WRITE(*,*) 'x(', i ,')*=',  x_solution(i)
                   END DO
               END IF
               WRITE(*,*) '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - '               
               WRITE(*,*) 'f*=', f_solution                    
               WRITE(*,*) '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - '   
               WRITE(*,*) 'The number of main iterations:', iter_counter
               WRITE(*,*) 'The number of subproblems solved:', subprob_counter             
               WRITE(*,*) 'The number of function values evaluated for DC component f_1:', f_counter 
               WRITE(*,*) 'The number of function values evaluated for DC component f_2:', f_counter 
               WRITE(*,*) 'The number of subgradients calculated for DC component f_1:', subgrad1_counter 
               WRITE(*,*) 'The number of subgradients calculated for DC component f_2:', subgrad2_counter 
               WRITE(*,*) 'the number of times approximate stopping condition was tested:', stop_cond_counter              
               WRITE(*,*) '---------------------------------------------------------------------------'     
               WRITE(*,*) 'CPU time used:', finish_time-start_time                          
               WRITE(*,*) '--------------------------------------------------------------------------- '
               WRITE(*,*) 'Elapsed time:', elapsed_time
               WRITE(*,*) '---------------------------------------------------------------------------'                
           END IF 

           CALL deallocate_b1(B1)       
           CALL deallocate_b2(B2) 
           
           END SUBROUTINE PBDC_algorithm      
        !.......................................................................................           
        ! <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>   
        ! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  
        !| | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |        
        !** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** 
        !|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|


        
        
        
        ! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  
        !| | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |
        ! START ** START ** START ** START ** START ** START ** START ** START ** START ** START 
        !|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|         
        ! <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>
        !***************************************************************************************
        !  ----------------------------------------------------------------------------------  |
        !  |                                                                                |  |
        !  |                             THE MAIN ITERATION                                 |  |
        !  |                                                                                |  |
        !  ----------------------------------------------------------------------------------  |
        !***************************************************************************************        
           
            SUBROUTINE PBDC_main_iteration( x_k, f1_k, f2_k, x_new, f1_new, f2_new,&
                            & f_0, B1, B2, L1, L2, crit_tol, eps, m, c, r_dec, r_inc,&
                            & mrounds, iprint, reason_for_stop, iter_counter, &
                            & subprob_counter, fi_counter, subgrad1_counter, subgrad2_counter,&
                            & stop_cond_counter, agg_in_use, stepsize_used, problem1, problem2, user_n)                 
            !
            ! Executes the 'main iteration' algorithm. It is needed to find the search direction at the current iteration point. 
            !
            ! INPUT: * 'x_k'               : The current iteratioin point
            !        * 'user_n'            : The dimension of the problem
            !        * 'problem1'          : The f1 of the problem
            !        * 'problem2'          : The f2 of the problem
            !        * 'f1_k' and 'f2_k'   : The values of DC components f_1 and f_2 at the current iteration point
            !        * 'f_0'               : The value of the objective function f at a starting point
            !
            !        * 'L1' and 'L2'       : Lipschitz constants of f_1 and f_2
            !        * 'crit_tol' and 'eps': The criticality tolerance and the proximity measure
            !        * 'm'                 : The descent parameter
            !        * 'c'                 : The decrease parameter         
            !        * 'r_dec' and 'r_inc' : The decrease and increase parameters   
            !
            !        * 'mrounds'           : The maximum number of possible rounds in the 'main iteration'
            !        * 'iprint'            : Specifies print
            !
            !        * 'agg_in_use'        : If .TRUE. then aggregation is used in the algorithm
            !        * 'stepsize_used'     : If .TRUE. then simple stepsize determination is used in at the end of each 'main iteration'
            !
            !
            ! OUTPUT: * 'x_new'              : the new iteration point
            !         * 'f1_new' and 'f2_new': the new value of DC components f_1 and f_2
            !         * 'reason_for_stop'    : Indicates the cause of the termination in the 'main iteration' 
            !
            !         * 'iter_counter'       : the number of rounds needed in 'main iteration'
            !         * 'subprob_counter'    : the number of subproblems solved in 'main iteration'
            !         * 'fi_counter'         : the number of function values evaluated for a DC component in 'main iteration' (same for f_1 and f_2)
            !         * 'subgrad1_counter'   : the number of subgradients calculated for f_1 in 'main iteration'
            !         * 'subgrad2_counter'   : the number of subgradients calculated for f_2 in 'main iteration'
            !         * 'stop_cond_counter'  : the number of times approximate stopping condition is tested during the 'main iteration' 
            !
            ! INOUT: * 'B1' and B2': the bundles of the DC components f_1 and f_2
            !
            ! NOTICE: The dimensions of vectors 'x_k' and 'x_new' has to be same. (i.e. the dimensio of 'x_k' and 'x_new' has to be
            !         'user_n' when SUBROUTINE PBDC_main_iteration() is used in SUBROUTINE bundle_method.
            !
            !***********************************************************************************
               IMPLICIT NONE
            !**************************** NEEDED FROM THE USER *********************************  
            
               TYPE(kimppu1), INTENT(INOUT) :: B1                   ! the bundle B_1 for the DC component f_1
               TYPE(kimppu2), INTENT(INOUT) :: B2                   ! the bundle B_2 for the DC component f_2
               REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: x_k      ! the current iteration point
               REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: x_new    ! the new iteration point if 'reason_for_stop'=0 in the 'main iteration'

               REAL(KIND=dp), INTENT(IN)  :: f1_k, f2_k      ! the value of f_1 and f_2 at x_k
               REAL(KIND=dp), INTENT(OUT) :: f1_new, f2_new  ! the value of f_1 and f_2 at x_new if 'reason_for_stop=0' in the 'main iteration'
               REAL(KIND=dp), INTENT(IN)  :: f_0             ! the value of f at the starting point x_0
               
               REAL(KIND=dp), INTENT(IN) :: L1, L2        ! Lipschitz constants of f_1 and f_2
               REAL(KIND=dp), INTENT(IN) :: crit_tol, eps ! crit_tol=criticality tolerance and eps=proximity measure (needed in stopping conditions)
               REAL(KIND=dp), INTENT(IN) :: m             ! descent parameter
               REAL(KIND=dp), INTENT(IN) :: c             ! c=decrease parameter               
               REAL(KIND=dp), INTENT(IN) :: r_dec, r_inc  ! r_dec=decrease parameter and R_inc=increase parameter
               
               INTEGER :: user_n                          ! the dimension
               INTEGER, DIMENSION(:) :: problem1          ! the f1
               INTEGER, DIMENSION(:) :: problem2          ! the f2
               
               INTEGER :: mrounds                         ! the maximum number of possible rounds in the 'main iteration'
                                                          ! If mrounds<=0, then DEFAULT value 500 is used (This is also done in the SUBROUTINE bundle_method()
               
               INTEGER, INTENT(OUT) :: reason_for_stop    ! 0 - a new iteration point found
                                                          ! 1 - stopping condition satisfied (criticality)
                                                          ! 2 - approximate stopping condition satisfied (eps-criticality)
                                                          ! 3 - maximum number of rounds executed in 'main iteration'
                                                          
               INTEGER, INTENT(OUT) :: iter_counter        ! the number of rounds used in 'main iteration'
               INTEGER, INTENT(OUT) :: subprob_counter     ! the number of subproblems solved in 'main iteration'
               INTEGER, INTENT(OUT) :: fi_counter          ! the number of function values evaluated for a DC component in 'main iteration' (same for f_1 and f_2)
               INTEGER, INTENT(OUT) :: subgrad1_counter    ! the number of subgradients calculated for f_1 in 'main iteration'
               INTEGER, INTENT(OUT) :: subgrad2_counter    ! the number of subgradients calculated for f_2 in 'main iteration'  
               INTEGER, INTENT(OUT) :: stop_cond_counter   ! the number of times approximate stopping condition is tested during the 'main iteration'  
               
               INTEGER, INTENT(IN) :: iprint    ! the variable that specifies print:
                                                !   iprint = 0 : print is suppressed
                                                !   iprint = 1 : basic print of final result (nothing is printed in main iteration)
                                                !   iprint = -1: basic print of final result without the solution vector (nothing is printed in main iteration)                                             
                                                !   iprint = 2 : extended print of final result (nothing is printed in main iteration)
                                                !   iprint = -2: extended print of final result without the solution vector (nothing is printed in main iteration)
                                                !   iprint = 3 : basic print of intermediate results and extended print of final results (nothing is printed in main iteration)
                                                !   iprint = -3: basic print of intermediate results and extended print of final results without the solution vector (nothing is printed in main iteration)
                                                !   iprint = 4 : extended print of intermediate and final results (something is printed about the 'main iteration' (this print is done in the bundle algorithm))               
                                                !   iprint = -4: extended print of intermediate and final results without the solution vector (something is printed about the 'main iteration' (this print is done in the bundle algorithm))               
                                                !   iprint = 5 : everything is printed step by step (this is the only option which causes print in the 'main iteration') 
       
               LOGICAL :: agg_in_use            ! If .TRUE. then aggregation is used in the algorithm   
               LOGICAL :: stepsize_used         ! If .TRUE. then simple stepsize determination is used at the end of the 'main iteration'
               

           !***************************** LOCAL VARIABLES ************************************
               
               REAL(KIND=dp), DIMENSION(give_n_b1(B1)) :: d_t                    ! the search direction
               REAL(KIND=dp), DIMENSION(give_n_b1(B1)) :: y                      ! the new auxilary point
               REAL(KIND=dp), DIMENSION(give_n_b1(B1)) :: new_grad1, new_grad2   ! the subgradients of f_1 and f_2 at y
               REAL(KIND=dp), DIMENSION(give_n_b1(B1)) :: old_grad               ! the overwritten subgradient of f_1
               REAL(KIND=dp), DIMENSION(give_n_b1(B1)) :: agg_grad               ! the new aggregated subgradient of f_1
               REAL(KIND=dp), DIMENSION(give_n_b1(B1)) :: apu_grad               ! 'help' subgradient
               REAL(KIND=dp), DIMENSION(give_n_B1(B1)) :: bxi                    ! the solution to the norm minimization problem (i.e. bxi = bxi_1^* - bxi_2^*)
               REAL(KIND=dp), DIMENSION(give_n_b1(B1)) :: vect                   ! 'help' vector               
               
               REAL(KIND=dp) :: norm1        ! ||\bxi_1(x_k)||
               REAL(KIND=dp) :: max_norm     ! the value of the maximum subgradient norm ||\bxi_{2,max}||
               REAL(KIND=dp) :: eps1         ! eps_1 needed in the 'main iteration'
               REAL(KIND=dp) :: t            ! the proximity parameter t
               REAL(KIND=dp) :: div_t        ! 1.0_dp divided by the proximity parameter t
               REAL(KIND=dp) :: t_min, t_max ! the bounds for proximity parameter t
               REAL(KIND=dp) :: theta        ! the proximity measure
               
               REAL(KIND=dp) :: d_norm             ! the norm of the direction vector d_t
               REAL(KIND=dp) :: delta1, delta2     ! the predicted changes of f_1 and f_2, respectively
               REAL(KIND=dp) :: f1_y, f2_y         ! the function values of f_1 and f_2 at y
               REAL(KIND=dp) :: f_k                ! the function value of f at x_k
               REAL(KIND=dp) :: real_decrease      ! the real decrease of the objective function f
               REAL(KIND=dp) :: lin_err1, lin_err2 ! the linearization errors of f_1 and f_2 (calculated using y)
               REAL(KIND=dp) :: old_lin_err        ! the overwritten linearization error of f_1
               REAL(KIND=dp) :: agg_lin_err        ! the new aggregated linearization error of f_1
               
               REAL(KIND=dp) :: obj              ! the objective function value of the norm minimization problem at the solution 'bxi'             
               REAL(KIND=dp) :: help             ! 'help' variable  

               INTEGER :: s_counter              ! the number of subproblems solved during the current iteration round of the 'main iteration' algorithm
               INTEGER :: N1_before, N2_before   ! the sizes of bundles before the removal of elements when the norm minimization problem is solved
               INTEGER :: N1, N2                 ! the sizes of bundles when we start to solve the norm minimization problem
               INTEGER :: i, j , ind 
               INTEGER :: it
                           
               LOGICAL :: t_changed              ! .TRUE. if t has been changed during the previous round (.TRUE. also when initialization is done)
               LOGICAL :: was_b1_full            ! .TRUE. if during the previous round something was added to B_1 and B_1 was full before this insertion
               LOGICAL :: agg_used               ! .TRUE. if aggregation element is in use
               LOGICAL :: stop_main_it           ! .TRUE. if the current 'main iteration' can be stopped
               
               LOGICAL :: stop_askelpituus       ! .TRUE. if the current 'main iteration' can be stopped
               REAL(KIND=dp) :: askelpituus        ! .TRUE. if the current 'main iteration' can be stopped
               REAL(KIND=dp) :: koe_f1            ! .TRUE. if the current 'main iteration' can be stopped
               REAL(KIND=dp) :: koe_f2            ! .TRUE. if the current 'main iteration' can be stopped
               REAL(KIND=dp), DIMENSION(give_n_b1(B1)) :: koe_y            ! .TRUE. if the current 'main iteration' can be stopped
                       

           !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^     
           ! <>  <>  <>  <>  <>  <>  <>  MAIN ITERATION STARTS  <>  <>  <>  <>  <>  <>  <>  <>
           !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
           
           !__________________________________________________________________________________     
           !************************ STEP 0: STOPPING CONDITION ******************************
                   
               iter_counter = 1                   ! the first round is executed
               subprob_counter = 0                ! the number of subproblems solved
               fi_counter = 0                     ! the number of function values evaluated
               subgrad1_counter = 0               ! the number of subgradients calculated for f_1
               subgrad2_counter = 0               ! the number of subgradients calculated for f_2
               stop_cond_counter = 0              ! the number of times approximate stopping condition is tested
               
               it = 0    
               agg_used = .FALSE.                 ! Aggregation is initialized to be .FALSE.
           
               stop_main_it = .FALSE.
               
               vect = give_subgrad_b1(B1,0) - give_subgrad_b2(B2,0)       ! vector: ' \bxi_1(x_k) - \bxi_2(x_k) '

               ! NOTICE: we use 0.1 * crit_tol              
               start: IF ( SQRT(DOT_PRODUCT(vect,vect)) < (0.1_dp*crit_tol)) THEN  ! ||\bxi_1(x_k) - \bxi_2(x_k)|| < 0.1 * crit_tol 
               !**>>**>>**>> APPROXIMATE CRITICALITY OBTAINED <<**<<**<<**
               
                   stop_main_it = .TRUE.  ! the 'main iteration' can be stopped
                   reason_for_stop = 1    ! the reason for stopping is criticality 

        
                   IF(iprint == 5) THEN  ! print everything        
                       WRITE(*,*) 'At STEP 0 the current point satisfies CRITICALITY'   
                   END IF

                   
           !__________________________________________________________________________________         
           !************************ STEP 0: END *********************************************         
               ELSE start   
               !_______________________________________________________________________________
               !************************ STEP 1: PARAMETER INITIALIZATION *********************
           
                   f_k = f1_k - f2_k                    ! the function value of f at x_k
                   
                   vect = give_subgrad_b1(B1,0)         ! the vector \bxi_1(x_k)
                   norm1 = SQRT(DOT_PRODUCT(vect,vect)) ! the value of the norm ||\bxi_1(x_k)||
                   max_norm = max_norm_value(B2)        ! the value of the maximum subgradient norm in the bundle B_2
                   eps1 = 0.5_dp * eps / ( MAX(L1,L2) ) ! the value of epsilon_1
                   t_min = (0.5_dp * r_dec * eps1) / ( norm1 + max_norm ) ! the lower bound for t
                   t_max = r_inc * t_min                                  ! the upper bound for t
                   theta = r_dec * t_min * crit_tol                       ! theta
                               
                   t = select_value_t(t_min,t_max)    ! the parameter t is selected         

                   ! print everything   
                   IF(iprint == 5) THEN  
                      WRITE(*,*) ' '                   
                      WRITE(*,*) 'STEP 1: INITIALIZATION OF PARAMETERS'
                      WRITE(*,*) ' '                   
                      WRITE(*,*) 't=', t
                      WRITE(*,*) 't is from the interval [t_min, t_max]=[', t_min, ',', t_max, ']' 
                      WRITE(*,*) 'theta=', theta, ' eps=', eps
                      WRITE(*,*) 'The maximum subgradient norm in B_2: ', max_norm, ' The norm |||\bxi_1(x_k)||: ', norm1 
                      WRITE(*,*) ' '
                      WRITE(*,*) '-------------------------------------------------------------------------------'  
                   END IF                  
            
            
                   t_changed = .TRUE.                 ! the parameter t was changed
                   was_b1_full = .FALSE.              ! B_1 was not full during the previous round since this is the first round of 'main iteration'
                       
                                   
                   IF ( mrounds <= 0 ) THEN           
                       mrounds = 5000                 ! the DEFAULT value for 'mrounds' is used
                   END IF
               !_______________________________________________________________________________
               !************************ STEP 1: END ****************************************** 
               
               END IF start

               
           !----------------------------------------------------------------------------------   
           ! <>  <>  <>  <>  <>  <>  <> REPEATABLE PART BEGINS <>  <>  <>  <>  <>  <>  <>  <>
           !----------------------------------------------------------------------------------
           
               DO WHILE ( (.NOT. stop_main_it ) .AND. (iter_counter <= mrounds ))   ! is repeated until the 'main iteration' can be terminated         
               !______________________________________________________________________________
               !************************ STEP 2: SEARCH DIRECTION ****************************
                   
               
                   IF (agg_in_use) THEN             ! 'agg_in_use' tells whether we use aggregation or not. If 'agg_in_use=.FALSE. then also 'agg_used'=.FALSE.
                       agg_used = is_agg_used(B1)   ! tells whether aggregation element is used or not? 
                   END IF   
                   
                   IF ( agg_used ) THEN    ! IF we have .FALSE. here then the algorithm does NOT use AGGREGATION

                        CALL subproblem_solver(x_k, give_n_b1(B1), give_size_b1(B1)+2, &  ! subproblems are solved with aggregation
                                          &  B1, B2, t, s_counter)
                   
                   ELSE
                   
                        CALL subproblem_solver(x_k, give_n_b1(B1), give_size_b1(B1)+1, &  ! subproblems are solved without aggregation
                                          &  B1, B2, t, s_counter)
                   END IF   
                    
                   subprob_counter = subprob_counter + s_counter    ! the number of subproblems solved so far
                                   
                   CALL add_glob_index(B2)              ! calculates the index of the subproblem yielding the global solution        
                   d_t = give_solution(B2)              ! the global solution d_t is selected   
                   d_norm = SQRT(DOT_PRODUCT(d_t,d_t))  ! the norm ||d_t|| of the solution d_t is calculated

               
                   ! Aggregated element is calculated. 
                   ! NOTICE: The aggregated element (if used) is now always added into bundle B_1 (also in those cases when the bundle B_1 is not full)

                   ! It is also possible to add the aggregated element into bundle B_1 only when the bundle B_1 is full
                   ! IF (is_full_b1(B1)) THEN

                   IF (agg_in_use) THEN                     ! Is done only when 'agg_in_use'=.TRUE.
                       ind = give_solution_ind(B2)
                       apu_grad = give_subgrad_b2(B2,ind)
                       div_t = 1.0_dp / t
                       DO i = 1, give_n_b2(B2)
                          agg_grad(i) = (-div_t) * d_t(i) + apu_grad(i)
                       END DO
                       agg_lin_err = - give_decrease(B2) - (div_t) * DOT_PRODUCT(d_t,d_t)
                       agg_lin_err = agg_lin_err + give_linerr_b2(B2,ind)
                       CALL add_agg_element_b1(B1,agg_grad,agg_lin_err)
                   END IF

                   !END IF  
           
                   ! print everything                   
                   IF(iprint == 5) THEN
                      WRITE(*,*) ' '
                      WRITE(*,*) 'STEP 2: SEARCH DIRECTION' 
                      WRITE(*,*) ' '                  
                      WRITE(*,*) 'Round:' , iter_counter  
                      WRITE(*,*) 'Search direction: [', d_t , ']' 
                      WRITE(*,*) 'The number of subproblems solved during this round: ', s_counter
                      WRITE(*,*) 'The number of subproblems solved till now: ', subprob_counter
                      WRITE(*,*) ' '                  
                      WRITE(*,*) '-----------------------------------------------------------------------------------'
                   END IF 
 
               !______________________________________________________________________________
               !************************ STEP 2: END *****************************************     
                   
               !->->->->->-> EXECUTION OF BRANCH BEGINS (2 POSSIBLE BRANCHES) <-<-<-<-<-<-<-<-
                   branches: IF (d_norm < theta) THEN               
               !->->->->->->->->->->->->->-> BRANCH 1 BEGIN <-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-             
                   !__________________________________________________________________________
                   !******************** STEP 3: APPROXIMATE STOPPING CONDITION **************
                       
                       N1_before = give_size_b1(B1) + 1
                       N2_before = give_size_b2(B2) + 1
                        
                       CALL delete_b1(B1, eps)    ! the bundle elements for which 'alpha_{1,j}>eps' are deleted from B_1
                       CALL delete_b2(B2, eps)    ! the bundle elements for which 'alpha_{2,j}>eps' are deleted from B_2

                       N1 = give_size_b1(B1) + 1  ! the current size of the bundle B_1 (NOTICE: the 'current bundle element' is also taken into account)
                       N2 = give_size_b2(B2) + 1  ! the current size of the bundle B_2 (NOTICE: the 'current bundle element' is also taken into account)
                       
                       
                       IF ( (N1_before == N1) .AND. (N2_before == N2)) THEN 
                          obj = d_norm / t
                       ELSE
                          ! The norm minimization problem is solved. The solution is bxi = bxi_1^* + bxi_2^* 
                          CALL norm_solver(bxi, obj, N1+N2, N1, N2, B1, B2, user_n)
                                           
                       END IF

                       stop_cond_counter = stop_cond_counter + 1   ! the update of the stopping condition counter
                       
                       IF ( obj < crit_tol ) THEN
                       !**>>**>>**>> APPROXIMATE EPS-CRITICALITY OBTAINED <<**<<**<<**
                           stop_main_it = .TRUE.  ! the 'main iteration' can be stopped
                           reason_for_stop = 2    ! the reason for stopping is eps-criticality
                           
                           ! print everything
                           IF(iprint == 5) THEN
                              WRITE(*,*) ' '
                              WRITE(*,*) 'STEP 3: APPROXIMATE STOPPING CONDITION'   
                              WRITE(*,*) ' '                  
                              WRITE(*,*) 'The number of the approximate stopping condition' , stop_cond_counter
                              WRITE(*,*) 'The solution of the norm problem "bxi"= [', bxi, ']' 
                              WRITE(*,*) 'The objective function of the norm problem:', obj                                    
                              WRITE(*,*) ' '
                              WRITE(*,*) 'The solution bxi satisfies the APPROXIMATE STOPPING CONDITION.'
                              WRITE(*,*) ' ' 
                           END IF 
                           
                       ELSE
                           t_max = t_max - r_dec * ( t_max - t_min ) ! the upper bound t_max is reduced
                           t = select_value_t(t_min,t_max)           ! the selection of the parameter t
                           t_changed = .TRUE.                        ! the parameter t was updated
                           iter_counter = iter_counter + 1           ! the update of the iteration counter  
                       
                           ! print everything
                           IF(iprint == 5) THEN
                              WRITE(*,*) ' '
                              WRITE(*,*) 'STEP 3: APPROXIMATE STOPPING CONDITION'   
                              WRITE(*,*) ' '                  
                              WRITE(*,*) 'The number of the approximate stopping condition' , stop_cond_counter
                              WRITE(*,*) 'The solution of the norm problem "bxi"= [', bxi, ']' 
                              WRITE(*,*) 'The objective function of the norm problem:', obj
                              WRITE(*,*) ' '                  
                              IF (obj >= crit_tol) THEN
                                 WRITE(*,*) 'A new upper bound t_max used and a new interval for t is: [', t_min, t_max, ']'  
                                 WRITE(*,*) 'Due to this a new parameter t=', t 
                              END IF
                              WRITE(*,*) ' '                  
                              WRITE(*,*) '-----------------------------------------------------------------------------------'
                           END IF                      
                       
                       END IF   
 

                   !__________________________________________________________________________
                   !******************** STEP 3: APPROXIMATE STOPPING CONDITION **************  
                
               !->->->->->->->->->->->->->-> BRANCH 1 END <-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-
                        
                   ELSE branches
               !->->->->->->->->->->->->->-> BRANCH 2 BEGIN <-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-
                   !__________________________________________________________________________
                   !******************** STEP 4: DESCENT TEST ********************************
       
                       y = x_k + d_t    ! a new auxilary point y
                       f1_y = f1(y,problem1,user_n)     ! the value of f_1 at y
                       f2_y = f2(y,problem2,user_n)     ! the value of f_2 at y
                       
                       fi_counter = fi_counter + 1 ! one more objective function value calculated for a DC component
                       
                       real_decrease = ( f1_y - f2_y ) -  f_k   ! the real decrease of f                   

                   ! print everything
                   IF(iprint == 5) THEN
                       WRITE(*,*) ' '
                       WRITE(*,*) 'STEP 4: DESCENT TEST'    
                       WRITE(*,*) ' '                 
                       WRITE(*,*) 'A new AUXILIARY point: [', y, ']'                  
                       WRITE(*,*) ' '                 
                       WRITE(*,*) 'Real decrease:' , real_decrease, ' Predicted decrease: ', give_decrease(B2)            
                       WRITE(*,*) ' Is DECREASE enough: ', real_decrease, '<=', (m * give_decrease(B2)), '?'                      
                       WRITE(*,*) ' ' 
                   END IF
                       
                   
                       branch2: IF ( real_decrease <= (m * give_decrease(B2)) ) THEN
                       !**>>**>>**>> NEW ITERATION POINT FOUND <<**<<**<<**  
                           
                           IF(stepsize_used) THEN 
                           !-**--**-  STEPSIZE DETERMINATION BEGIN -**--**-
                           IF (real_decrease < -0.1_dp) THEN
                              askelpituus = 1.0_dp
                              stop_askelpituus = .FALSE.
                              DO WHILE((.NOT. stop_askelpituus))
                                 koe_y = y + d_t
                                 koe_f1 = f1(koe_y,problem1,user_n)     ! the value of f_1 at koe_y
                                 koe_f2 = f2(koe_y,problem2,user_n)     ! the value of f_2 at koe_y 
                                 fi_counter = fi_counter + 1    
                                 IF ( (koe_f1 - koe_f2 - f_k) <= ( 0.7_dp * m * give_decrease(B2))) THEN 
                                    askelpituus = askelpituus + 1.0_dp
                                    y = koe_y
                                    f1_y = koe_f1
                                    f2_y = koe_f2                                   
                                 ELSE
                                    stop_askelpituus = .TRUE.  
                                 END IF
                              END DO
                           END IF 
                           !-**--**- STEPSIZE DETERMINATION END -**--**-
                           END IF
                           
                           x_new = y               ! the new iteration point is the current auxilary point
                           f1_new = f1_y           ! the value of f_1 at x_new
                           f2_new = f2_y           ! the value of f_2 at x_new
                           stop_main_it = .TRUE.   ! the 'main iteration' can be stopped
                           reason_for_stop = 0     ! the reason for stopping is the new iteration point
                           
                           ! print everything
                           IF(iprint == 5) THEN
                              WRITE(*,*) 'A NEW ITERATION POINT is found: x_new= [', x_new, ']'  
                              WRITE(*,*) 'A NEW OBJECTIVE VALUE: f_new=', f1_new - f2_new                 
                              WRITE(*,*) 'The old point: [', x_k,']'              
                              WRITE(*,*) 'The old objective value:', f_k              
                              WRITE(*,*) ' '                  
                           END IF
                           
                           
                   !__________________________________________________________________________
                   !******************** STEP 4: END *****************************************         
                   
                       ELSE branch2 
                       !----------- SOME UPDATES IN VARIABLES BEGINS -------------------------  

                          ! print everything                       
                          IF (iprint == 5) THEN
                              WRITE(*,*) '-----------------------------------------------------------------------------------'
                          END IF                       
                           
                       
                   !__________________________________________________________________________     
                   !******************** STEP 5: BUNDLE UPDATE *******************************
                          
                          update: IF ( ((f1_y - f2_y - f_0)>0 ) .AND. & 
                                                      & (d_norm > eps1)) THEN 
                          !----------- REDUCTION OF PARAMETER t --------------------------------
                              t = t - r_dec * ( t - t_min )  ! the parameter t is updated
                              t_changed = .TRUE.             ! the parameter t was changed
                

                          ! print everything
                          IF(iprint == 5) THEN
                              WRITE(*,*) ' '
                              WRITE(*,*) 'STEP 5a: PARAMETER t UPDATE'  
                              WRITE(*,*) ' '                  
                              WRITE(*,*) 'A new parameter t=', t
                              WRITE(*,*) 'Reason for this is that conditions:'
                              WRITE(*,*)  '(f1_y - f2_y - f_0)>0 ) .AND. (d_norm > eps1) are fulfilled:'
                              WRITE(*,*) ' '                  
                              WRITE(*,*) '-----------------------------------------------------------------------------------'      
                          END IF                  

                               
                          ELSE update
                          !----------- BUNDLE AND PARAMETER UPDATE --------------------------
                              
                              i = give_solution_ind(B2)       ! the index of the subproblem yielding the global solution
                              vect = give_subgrad_b2(B2, i)   ! the subgradient of f_2 in the subproblem yielding the global solution
                              delta2 = - DOT_PRODUCT(vect, d_t) + give_linerr_b2(B2,i)  ! the value of delta_2
                              delta1 = give_decrease(B2) - delta2                       ! the value of delta_1
                              
                              new_grad1 = subgradient_f1(y,problem1,user_n )                        ! a subgradient of f_1 at y
                              subgrad1_counter = subgrad1_counter + 1              ! a new subgradient was evaluated for f_1            

                              lin_err1 = f1_k - f1_y + DOT_PRODUCT(d_t,new_grad1)  ! a new linearization error for f_1
                              
                              IF (real_decrease > -m * give_decrease(B2) .OR. it > 20)  THEN    ! adjustment of the proximity parameter if necessary
                                  t =  t - c * ( t - t_min )
                              END IF          
                             
                            !->>>>>>>>>>>> BUNDLE B_1 UPDATE BEGINS <<<<<<<<<<<<<<<<<<<<<<<<<<- 
                              bundle1: IF (is_full_b1(B1)) THEN   
                              ! bundle B_1 is full and a new element is added to the bundle B_1
                              ! the overwritten element is written/taken down 
                                  CALL add_element_b1(B1, new_grad1, lin_err1)
                                  !CALL add_give_b1(B1, new_grad1, lin_err1, old_grad,old_lin_err)
                                  was_b1_full = .TRUE.
                                
                              ELSE bundle1
                              ! bundle B_1 is NOT full and a new element is added to the bundle B_1
                                  CALL add_element_b1(B1, new_grad1, lin_err1)
                                  was_b1_full = .FALSE.

                              END IF bundle1                             
                            !->>>>>>>>>>>> BUNDLE B_1 UPDATE ENDS <<<<<<<<<<<<<<<<<<<<<<<<<<<-
                            
                              ! print everything 
                              IF(iprint == 5) THEN
                                  WRITE(*,*) ' '
                                  WRITE(*,*) 'STEP 5b: BUNDLE UPDATE'   
                                  WRITE(*,*) ' '                  
                                  WRITE(*,*) 'delta2: ', delta2, ' delta1: ', delta1, ' glob solution index in B2: ', i 
                                  WRITE(*,*) ' ' 
                                  WRITE(*,*) 'A new bundle element is inserted into the bundle B_1. '   
                              END IF
  
                        
                            !->>>>>>>>>>>> BUNDLE B_2 AND PARAMETER UPDATE BEGINS <<<<<<<<<<<-             

                            bundle2: IF( delta2 >= 0 .AND. (give_max_size_b2(B2)>0) ) THEN   ! bundle B_2 is updated (happens only when delta2 >= 0 and bundle size is larger than 1) (SOME OTHER INSERTION RULES ALSO POSSIBLE!)
                                                             
                             
                                  new_grad2 = subgradient_f2(y,problem2,user_n)               ! a subgradient of f_2 at y  
                                  subgrad2_counter = subgrad2_counter + 1              ! a new subgradient was evaluated for f_2
                                  lin_err2 = f2_k - f2_y + DOT_PRODUCT(d_t,new_grad2)  ! a new linearization error for f_2  
                                  CALL add_element_b2(B2, new_grad2, lin_err2)         ! a new element is inserted into the bundle B_2
                                  !In the algorithm we never overwrite the 'bundle element' yielding the previous global solution of the search direction problem.

                              
                                  ! print everything                                  
                                  IF(iprint == 5) THEN
                                  WRITE(*,*) 'A new bundle element is inserted into the bundle B_2'
                                  WRITE(*,*) ' '                  
                                  WRITE(*,*) '-----------------------------------------------------------------------------------'
                                  END IF 

                                  
                   !__________________________________________________________________________     
                   !******************** STEP 5: END *****************************************
                   
                   !__________________________________________________________________________     
                   !******************** STEP 6: PARAMETER UPDATE ****************************
                   
                                  help = SQRT(DOT_PRODUCT(new_grad2,new_grad2)) ! the norm of the new subgradient of f_2 (subgradient is calculated 
                                                                                ! at the new auxilary point y)
                                  IF(help > max_norm) THEN 
                                      max_norm = help            ! the updated value of the maximum subgradient norm in the bundle B_2
                                      t_min = (0.5_dp * r_dec * eps1) / ( norm1 + max_norm ) ! the updated value of the lower bound t_min
                                      theta = r_dec * t_min * crit_tol                       ! the updated value of the parameter theta
                                  END IF            
                                  
                                  ! print everything
                                  IF(iprint == 5) THEN
                                  WRITE(*,*) ' '
                                  WRITE(*,*) 'STEP 6: PARAMETER UPDATE'           
                                  WRITE(*,*) ' '                  
                                  WRITE(*,*) '-----------------------------------------------------------------------------------'
                                  END IF 
  
                                  
                   !__________________________________________________________________________     
                   !******************** STEP 6: END *****************************************
                   
                            END IF bundle2


                            !->>>>>>>>>>>> BUNDLE B_2 AND PARAMETER UPDATE ENDS <<<<<<<<<<<<<-
                               
                              t_changed = .FALSE.   ! the parameter t was NOT changed during this ELSE branch of 'update'
                                  
                          END IF update        
                        !------ PARAMETER t REDUCTION & BUNDLE AND PARAMETER UPDATE ENDS -----
                        
                          iter_counter = iter_counter + 1  ! update of the iteration counter                        
                              
                       END IF branch2   
                     !--------- NEW ITERATION POINT & SOME UPDATES IN VARIABLES ENDS --------- 
                     
               !->->->->->->->->->->->->->-> BRANCH 2 END <-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-   
                   END IF branches
               !->->->->->->-> EXECUTION OF BRANCH ENDS  (2 BRANCHES) <-<-<-<-<-<-<-<-<-<-<-<-
              END DO 
              
           !----------------------------------------------------------------------------------   
           ! <>  <>  <>  <>  <>  <>  <>  REPEATABLE PART ENDS  <>  <>  <>  <>  <>  <>  <>  <> 
           !----------------------------------------------------------------------------------           
              
              IF ( (.NOT. stop_main_it ) ) THEN 
                  reason_for_stop = 3            ! the maximum number of rounds have been executed 
              END IF        
              
           END SUBROUTINE PBDC_main_iteration
        !.......................................................................................           
        ! <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>   
        ! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  
        !| | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |        
        !** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** 
        !|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|
        
            
        
        
        
        ! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  
        !| | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |
        ! START ** START ** START ** START ** START ** START ** START ** START ** START ** START  
        !|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|        
        ! <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>       
        !***************************************************************************************
        !  ----------------------------------------------------------------------------------  |
        !  |                                                                                |  |        
        !  |                           PARAMETER t SELECTION                                |  |
        !  |                                                                                |  |        
        !  ----------------------------------------------------------------------------------  |
        !***************************************************************************************           
           
           FUNCTION select_value_t(t_min,t_max) RESULT(t)
                IMPLICIT NONE
                REAL(KIND=dp), INTENT(IN) :: t_min, t_max
                REAL(KIND=dp) :: t, a
            
                t = 0.8_dp * (t_min + t_max) ! selection of t
                
                !OTHER OPTIONS:
                
!               t = 0.5_dp * (t_min + t_max)    ! selects the middle point from the interval [t_min, t_max]
!               t = t_min + 0.8*(t_max - t_min)  

!               t = t_min                   ! selects the lower bound
!               t = t_max                   ! selects the upper bound
                
           END FUNCTION select_value_t
        !.......................................................................................   
        ! <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>
        ! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  
        !| | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |        
        !** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** END **  
        !|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|    


  


    
          
        ! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  
        !| | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |
        ! START ** START ** START ** START ** START ** START ** START ** START ** START ** START 
        !|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|         
        ! <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>
        !***************************************************************************************
        !  ----------------------------------------------------------------------------------  |
        !  |                                                                                |  |        
        !  |                           SUBPROBLEM SOLVER                                    |  |
        !  |                                                                                |  |        
        !  ----------------------------------------------------------------------------------  |
        !***************************************************************************************           

            SUBROUTINE subproblem_solver(x, NF, NA, B1, B2, t,subprob_counter) 
               ! Solves the subproblems in the original search direction problem
               ! Calls PLQDF1 by Ladislav Luksan
               IMPLICIT NONE               
            !**************************** NEEDED FROM USER ***********************************             
 
               TYPE(kimppu1), INTENT(IN)    :: B1                  ! bundle B_1
               TYPE(kimppu2), INTENT(INOUT) :: B2                  ! bundle B_2
               REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: x        ! current iteration point
               
               REAL(KIND=dp), INTENT(IN) :: t                      ! proximity parameter t
               
               INTEGER, INTENT(IN)  :: NF    ! number of variables
               INTEGER, INTENT(IN)  :: NA    ! bundle size of B1         IF ( NA = give_size_b1(B1) + 2 ) THEN aggregation is used
               
               INTEGER, INTENT(OUT) :: subprob_counter    ! number of subproblems solved    
                       
              
            !*************************** LOCAL VARIABLES ************************************
            
               ! .. Parameters used in PLQDF1 ..
               REAL(KIND=dp), PARAMETER :: ETA0  = 1.0E-15_dp  
               REAL(KIND=dp), PARAMETER :: ETA2  = 1.0E-12_dp
               REAL(KIND=dp), PARAMETER :: ETA9  = 1.0E+60_dp
               REAL(KIND=dp), PARAMETER :: EPS7  = 1.0E-14_dp
               REAL(KIND=dp), PARAMETER :: EPS9  = 1.0E-12_dp 
               
               ! .. Scalar Arguments in PLQDF1 ..               
               REAL(KIND=dp) :: GMAX   ! output of PLQDF1: maximum absolute value of a partial derivative
               REAL(KIND=dp) :: UMAX   ! output of PLQDF1: maximum absolute value of a negative lagrange multiplier
               REAL(KIND=dp) :: XNORM  ! output of PLQDF1: value of linearized minimax function = delta_1 + delta_2  
               
               REAL(KIND=dp) :: alpha_b2     ! a linearization error of B_2
               REAL(KIND=dp) :: alpha_b1     ! a linearization error of B_1
               REAL(KIND=dp) :: u, obj       ! 'help' variable
               
               REAL(KIND=dp) :: a
               INTEGER :: i, j, k, l
               
               INTEGER, PARAMETER :: NC = 0             ! number of constraints is zero
               INTEGER :: IDECF=10, KBC=0, KBF=0, MFP=2 ! IDECF=10 diagonal matrix; KBC=0 no linear constraints; KBF=0 no simple bounds; MFP=2 optimum feasible point
               INTEGER :: ITERQ, N                      ! output values of PLQDF1: ITERQ=type of feasible point 
                                                        ! N=dimension of manifold defined by active constraints 
               ! .. Array Arguments in PLQDF1..
               INTEGER, DIMENSION(NF) :: IX     ! vector containing types of bounds
               INTEGER, DIMENSION(NA) :: IA     ! vector containing types of deviations 
               INTEGER, DIMENSION(NC) :: IC     ! vector containing types of constraints. NOT significant because NC=0
               INTEGER, DIMENSION(NF+1) :: IAA  ! Output of PLQDF1: vector containing indicies of active functions
               
               REAL(KIND=dp), DIMENSION(NF) :: XL, XU          ! lower and upper bounds for x. NOT significant because variables are unbounded
               REAL(KIND=dp), DIMENSION(NF) :: H               ! diagonal matrix (1/t)*I
               REAL(KIND=dp), DIMENSION(NC) :: CF, CL, CU      ! NOT significant since NC=0 (are related to constraints)
               REAL(KIND=dp), DIMENSION(NA) ::  AF             ! vector of bundle function values (-alpha)
               REAL(KIND=dp), DIMENSION(NA) :: AFD             ! Output of PLQDF1: vector containing increments of the approximated functions
               REAL(KIND=dp), DIMENSION(NF+1) :: AZ            ! Output of PLQDF1: vector of Lagrange multipliers
               REAL(KIND=dp), DIMENSION(NF+1) :: S             ! Output of PLQDF1: direction vector
               REAL(KIND=dp), DIMENSION(NF) :: direction       ! Actual direction vector used (notice dimension)
               REAL(KIND=dp), DIMENSION(NF+1) :: G             ! Output of PLQDF1: gradient of the Lagrangian function
               REAL(KIND=dp), DIMENSION(NF*NA) :: AG           ! matrix whose columns are bundle subgradients 
               REAL(KIND=dp), DIMENSION((NF+1)*(NF+2)/2) :: AR ! Output of PLQDF1: triangular decomposition of kernel of the orthogonal projection
               REAL(KIND=dp), DIMENSION(NF*NC) :: CG           ! NOT significant since NC=0. matrix whose columns are normals of the linear constraints            
               
               ! .. Some other varibles ..    
               REAL(KIND=dp), DIMENSION(NF*NA) :: grad_m_b1           ! subgradient matrix of B_1
               REAL(KIND=dp), DIMENSION(NA) :: alpha_m_b1             ! linearization error matrix of B_1
               REAL(KIND=dp), DIMENSION(NF) :: grad_b2                ! a subgradient of B_2
               REAL(KIND=dp), DIMENSION(NF) :: grad_b1                ! a subgradient of B_1


           !************************** SUBPROBLEM SOLVER STARTS *********************************   
           
           !****************************** INITIALIZATIONS **************************************
           
               IX = 0             ! types of bounds: 0 - unbounded variables
               IA = 2             ! types of deviations
               
               u = 1.0_dp / t
               H = u              ! diagonal matrix
               
               IF ( (give_size_b1(B1)+2) == NA) THEN     ! if this is TRUE then aggregation is in use
                   grad_m_b1 = grad_matrix_agg(B1)       ! subgradient matrix of B_1 with aggregation
                   alpha_m_b1 = lin_error_matrix_agg(B1) ! linearization error matrix of B_1 with aggregation
               ELSE
                   grad_m_b1 = grad_matrix(B1)       ! subgradient matrix of B_1
                   alpha_m_b1 = lin_error_matrix(B1) ! linearization error matrix of B_1
               END IF
               
               subprob_counter = give_size_b2(B2)   

               !$OMP PARALLEL DO PRIVATE(grad_b2,alpha_b2,direction,a,obj) & 
               !$OMP FIRSTPRIVATE(NA,X,IX,XL,XU,AFD,IA,IAA) &
               !$OMP FIRSTPRIVATE(AR,AZ,CF,IC,CL,CU,CG,G,H,S,KBF) &
               !$OMP FIRSTPRIVATE(KBC,XNORM,UMAX,GMAX,N,ITERQ) &
               !$OMP PRIVATE(i, AG, AF ) &
               !$OMP SHARED(grad_m_b1,alpha_m_b1,B2)               
                    
               !->->->->->->->->->->->->->-> EACH SUBPROBLEM IS SOLVED: BEGIN <-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-
                    subproblems1: DO i = 0, subprob_counter   ! each subproblem is looked through
                                                          
                        grad_b2  = give_subgrad_b2(B2, i) ! subgradient of B_2 in the subproblem i
                        alpha_b2 = give_linerr_b2(B2, i) ! linearization error of B_2 in the subproblem i   
                        
                        !$OMP CRITICAL  
                        DO j = 1, NA
                           k = (j-1)*NF
                           DO l = 1, NF
                               AG(k+l) = grad_m_b1(k+l) - grad_b2(l)
                           END DO
                           AF(j) = - alpha_m_b1(j) + alpha_b2 
                        END DO
                        !$OMP END CRITICAL  
                        
                        !Calls PLQDF1 by Ladislav Luksan
                        CALL PLQDF1(NF,NA,NC,X,IX,XL,XU,AF,AFD,IA,IAA, &
                              & AG,AR,AZ,CF,IC,CL,CU,CG,G,H,S,2,KBF, &
                              & KBC,10,ETA0,ETA2,ETA9,EPS7,EPS9,XNORM, &                        
                              & UMAX,GMAX,N,ITERQ)
    
                        DO j = 1, NF
                           direction(j) = S(j) 
                        END DO
              
                        a = DOT_PRODUCT(direction, direction)

                        obj =   XNORM + (a * u) / 2
                        !$OMP CRITICAL
                        CALL add_solution(B2, i , direction, XNORM, obj )   
                        !$OMP END CRITICAL
                    
                    END DO subproblems1
               !->->->->->->->->->->->->->-> EACH SUBPROBLEM IS SOLVED: END <-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-                      
               !$OPM END PARALLEL DO

               subprob_counter = subprob_counter + 1
               
               
           END SUBROUTINE subproblem_solver
        !.......................................................................................
        ! <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>        
        ! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  
        !| | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |        
        !** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** 
        !|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|







      END MODULE pbdc