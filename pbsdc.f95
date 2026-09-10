      MODULE pbsdc
      
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !| |                                                                                  | |
        !| |    PBSDC - THE PROXIMAL BUNDLE METHOD FOR SEPARABLE NONSMOOTH DC OPTIMIZATION    | | 
        !| |                                 (version 1)                                      | |
        !| |                                                                                  | |
        !| |                      by Kaisa Joki and Jenni Lampainen                           | |
        !| |                                                                                  | |
        !| |      Features :                                                                  | |
        !| |                                                                                  | |
        !| |           * Possibility to use either approximate criticality or                 | |
        !| |             stationarity in the stopping condition.                              | |
        !| |                                                                                  | |
        !| |           * During each round of 'main iteration' it is possible to              | |
        !| |             utilize OpenMP to calculate subproblems in parallel.                 | |
        !| |             In order to do this, you need to use '-fopenmp' in Makefile.         | |
        !| |                                                                                  | |                       
        !| |                                                                                  | |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*     
        !|                                                                                      |
        !|    The code utilizes new version of PLQDF1 by Ladislav Luksan as a quadratic solver. |
        !|                                                                                      |
        !|    The code utilizes PVMM by Ladislav Luksan as a norm minimization solver. This     | 
        !|    subroutine uses PQSUBS and MQSUBS by Ladislav Luksan. PVMM is a VARIABLE METRIC   |
        !|    ALGORITHM for UNCONSTRAINED and LINEARLY CONSTRAINED OPTIMIZATION.                |
        !|                                                                                      |
        !|    The subroutine PVMM together with PQSUBS and MQSUBS is licensed by                |
        !|    the GNU Lesser General Public License (LGPL).                                     |
        !|                                                                                      |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        
        !$USE omp_lib
        
        USE constants, ONLY : dp    ! Double precision (i.e. accuracy)
        USE bundle1                 ! The BUNDLE of the DC component f_1
        USE bundle2                 ! The BUNDLE of the DC component f_2
        USE functions               ! Contains INFORMATION from the USER
        USE norm_min                ! The SOLVER for the norm minimization problem
        
        IMPLICIT NONE    
        
        EXTERNAL PLQDF1             ! The QUADRATIC SOLVER by Ladislav Luksan
        
        
        TYPE block
           PRIVATE
           REAL(KIND=dp), DIMENSION(:), POINTER :: vector
        END TYPE block
        
        
        CONTAINS
                
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !| |                                                                                  | |
        !| |                            CONTAINS SUBROUTINES:                                 | | 
        !| |                                                                                  | |
        !| |      PBSDC_algorithm            : The proximal bundle method for block-separable | |
        !| |                                   nonsmooth DC optimization (i.e. PBSDC).        | | 
        !| |      PBSDC_main_iteration       : The main iteration algorithm needed in the     | |
        !| |                                   PBSDC.                                         | |
        !| |      PBSDC_guaranteeing_clarke  : The algorithm guaranteeing Clarke stationarity | |       
        !| |                                   for a solution obtained.                       | |       
        !| |      PBSDC_quadratic_solver     : The solver for the quadratic norm minimization | |       
        !| |                                   problem. (Needed in the Clarke                 | |
        !| |                                   stationary algorithm)                          | |       
        !| |      PBSDC_subproblem_solver    : The solver for subproblems in the search       | |
        !| |                                   direction problem. (Needed in the main         | |
        !| |                                   iteration algorithm)                           | |   
        !| |                                                                                  | |   
        !| |                                                                                  | |       
        !| |                            CONTAINS FUNCTIONS:                                   | |
        !| |                                                                                  | |
        !| |      PBSDC_select_value_t : Selects the parameter t from interval [t_min,t_max]. | |               
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
        !  | PBSDC - THE PROXIMAL BUNDLE ALGORITHM FOR SEPARABLE NONSMOOTH DC OPTIMIZATION  |  |
        !  |                                                                                |  |
        !  ----------------------------------------------------------------------------------  |
        !***************************************************************************************
        
           SUBROUTINE PBSDC_algorithm( f_solution, f0, mit,&
                            & mrounds, mrounds_clarke, termination, counter, CPUtime, agg_used,  &
                            & stepsize_used, iprint, problem10, problem20, user_n, n_blocks, size_blocks0, &
                            & starting_points, start_ind, time_limit, optimality_condition, boost, &
                            & nochange_mainit, nochange_tol )
                            
            ! 
            ! Solves the unconstrained nonsmooth DC minimization problem
            !
            ! INPUT:  * 'problem10'      : The f1 in problem
            !         * 'problem20'      : The f2 in problem     
            !         * 'user_n'         : The dimension of the problem    
            !         * 'n_blocks'       : The number of blocks   
            !         * 'size_blocks0'   : The sizes of the blocks   
            !         * 'mit'            : The maximum number of 'main iterations'      
            !         * 'mrounds'        : The maximum number of rounds during one 'main iteration'
            !         * 'mrounds_clarke' : The maximum number of rounds during one 'Clarke stationary' algorithm
            !         * 'agg_used'       : If .TRUE. then aggregation is used in the algorithm          
            !         * 'stepsize_used'  : If .TRUE. then simple stepsize determination is used in the algorithm            
            !         * 'iprint'         : Specifies the print          
            !         * 'starting_points': Table containing different starting points          
            !         * 'start_ind'      : The index of the used starting point          
            !         * 'time_limit'     : The time limit for the execution of the algorithm          
            !         * 'optimality_condition'  : The optimality condition used          
            !         * 'boost'          : 1 = no boosting, 2 = boosting is used       
            !         * 'nochange_mainit': 0=not used, integer value 'k' = algorithm is stopped if after 'k' main iterations change in the objective value is less than 'nochange_tol'      
            !         * 'nochange_tol'   : tolerance for change in the objective value          
            !
            ! OUTPUT: * 'f_solution'  : The objective function value at the solution 'x_solution' ("x_solution" is not returned for the user)
            !         * 'f0'          : The objective function value at the starting point
            !         * 'termination' : The cause of the termination in the algorithm
            !         * 'counter'     : Gives the values of different counters
            !         * 'CPUtime'     : the CPU time
            !
            ! NOTICE: * The dimension of the vector 'counter' has to be 8.
            !         * 'mit', 'mrounds' and 'mrounds_clarke' have to be integers.
            !         * IF ('mit' <= 0) THEN DEFAULT value 1000 is used
            !         * IF ('mrounds' <= 0 ) THEN DEFAULT value 5000 is used.
            !         * IF ('mrounds_clarke' <= 0 ) THEN DEFAULT value 5000 is used.
            !         * 'iprint' has to be -4, -3, -2, -1, 0, 1, 2, 3 or 4 . If it is NOT then DEFAULT value 1 is used. 
            !
            !***********************************************************************************
               IMPLICIT NONE
            !**************************** NEEDED FROM USER ************************************* 
              
               REAL(KIND=dp), DIMENSION(:,:), INTENT(IN) :: starting_points  ! the table containing different starting points
               INTEGER, INTENT(IN)                       :: start_ind        ! the index of the used starting point
               INTEGER, INTENT(INOUT)                    :: boost            ! 1 = no boosting, 2 = boosting is used
               
               REAL(KIND=dp), INTENT(OUT) :: f_solution    ! the objective function value at the solution 'x_solution'
               REAL(KIND=dp), INTENT(OUT) :: f0            ! the objective function value at the starting point 
               REAL(KIND=dp), INTENT(OUT) :: CPUtime       ! the CPU time
               REAL(KIND=dp), INTENT(IN) :: time_limit     ! the time limit for the execution of the algorithm 
               
               INTEGER, INTENT(IN) :: optimality_condition               ! the used optimality condition
               INTEGER, DIMENSION(n_blocks), INTENT(IN) :: problem10     ! the f1 in problem
               INTEGER, DIMENSION(n_blocks), INTENT(IN) :: problem20     ! the f2 in problem
               INTEGER, INTENT(IN) :: user_n                             ! the dimension of the problem
               INTEGER, INTENT(IN) :: n_blocks                           ! the number of blocks
               INTEGER, DIMENSION(n_blocks), INTENT(IN) :: size_blocks0  ! the vector containing the sizes of the blocks
               
               INTEGER, INTENT(IN) :: nochange_mainit       ! 0=not used, integer value 'k' = algorithm is stopped if after 'k' main iterations change in the objective value is less than 'nochange_tol'
               REAL(KIND=dp), INTENT(IN) :: nochange_tol    ! tolerance for change in the objective value
               
               INTEGER, INTENT(INOUT) :: mit               ! the maximum number of 'main iterations'
               INTEGER, INTENT(INOUT) :: mrounds           ! the maximum number of rounds during one 'main iteration'
               INTEGER, INTENT(INOUT) :: mrounds_clarke    ! the maximum number of rounds during one 'Clarke stationary' algorithm
                                                          
               INTEGER, INTENT(OUT) :: termination      ! Tells the reason for termination   
                                                        !-------------------------------------------------------------------------------------------------     
                                                          ! When stationarity condition is used
                                                             ! 1 - the stopping condition is satisfied (i.e. Clarke stationarity)
                                                             ! 2 - the approximate stopping condition is satisfied (i.e. the step-length beta* < eps)
                                                             ! 3 - the maximum number 'mrounds' of rounds executed in one main iteration
                                                             ! 4 - the maximum number of 'main iterations' is executed  
                                                             ! 5 - the maximum number 'mrounds_clarke' of rounds executed in one 'Clarke stationary' alqorithm
                                                             ! 6 - the time limit is met during the execution of the alqorithm
                                                             ! 7 - after 'nochange_mainit' main iterations change in the objective value is less than 'nochange_tol'
                                                        !-------------------------------------------------------------------------------------------------      
                                                          ! When criticality condition is used
                                                             ! 1 - the stopping condition is satisfied (i.e. criticality)
                                                             ! 2 - the approximate stopping condition is satisfied (i.e. eps-criticality)
                                                             ! 3 - the maximum number 'mrounds' of rounds is executed in one main iteration
                                                             ! 4 - the maximum number of 'main iterations' is executed 
                                                             ! 5 - the execution time exceeded the time limit
                                                             ! 7 - after 'nochange_mainit' main iterations change in the objective value is less than 'nochange_tol'
                                                         
               INTEGER, DIMENSION(8), INTENT(OUT) :: counter  ! contains the values of different counteres: 
                                                              !--------------------------------------------------------------------------------------------------------------------------       
                                                              ! When stationarity condition is used
                                                              !   counter(1) = iter_counter         the number of 'main iterations' executed
                                                              !   counter(2) = subprob_counter      the number of subproblems solved
                                                              !   counter(3) = f_counter            the number of function values evaluated for DC component in 'main iteration'
                                                              !   counter(4) = subgrad1_counter     the number of subgradients calculated for f_1 in 'main iteration'
                                                              !   counter(5) = subgrad2_counter     the number of subgradients calculated for f_2 in 'main iteration'    
                                                              !   counter(6) = stop_cond_counter    the number of times 'Clarke stationary algorithm' is executed 
                                                              !   counter(7) = clarke_f_counter     the number of function values evaluated for f in 'Clarke stationary algorithms'
                                                              !   counter(8) = clarke_sub_counter   the number of subgradients caluculated for f in 'Clarke stationary algorithms'
                                                              !--------------------------------------------------------------------------------------------------------------------------       
                                                              ! When criticality condition is used
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
           
           
           !***************************** LOCAL VARIABLES ************************************  
           
               TYPE(kimppu1), DIMENSION(:), ALLOCATABLE :: B1_table    ! The bundles B_1 for the block-saparable DC components f_1
               TYPE(kimppu2), DIMENSION(:), ALLOCATABLE :: B2_table    ! The bundles B_2 for the block-saparable DC components f_2

               REAL(KIND=dp), DIMENSION(user_n) :: x_start       ! the starting point (the dimension 'user_n' is the number of variables)
               REAL(KIND=dp), DIMENSION(user_n) :: x_solution    ! the solution obtained to the problem
               
               REAL(KIND=dp), DIMENSION(user_n) :: vect          ! help vector
               
               REAL(KIND=dp), DIMENSION(user_n) :: x_current     ! the current iteration point (the dimension 'user_n' is the number of variables)
               REAL(KIND=dp), DIMENSION(user_n) :: x_new         ! the new iteration point (obtained from the previous 'main iteration')              
                 
               TYPE(block), DIMENSION(:), ALLOCATABLE :: grad1, grad2  ! subgradients (the dimension 'user_n' is the length of subgradient)       
               TYPE(block), DIMENSION(:), ALLOCATABLE :: dd            ! the search direction obtained from the 'main iteration': dd=x_new -x_current 

               INTEGER, DIMENSION(n_blocks) :: problem1          ! the f1 in problem
               INTEGER, DIMENSION(n_blocks) :: problem2          ! the f2 in problem
               INTEGER, DIMENSION(n_blocks) :: size_blocks       ! the vector containing the sizes of the blocks
               REAL(KIND=dp), DIMENSION(user_n) :: x_0           ! the starting point (the dimension 'user_n' is the number of variables)           
     
               REAL(KIND=dp) :: L1, L2         ! Lipschitz constants of DC components f_1 and f_2

               REAL(KIND=dp) :: crit_tol, eps  ! 'crit_tol'=the stopping tolerance and 'eps'=the enlargement parameter
               REAL(KIND=dp) :: m              ! the descent parameter used in 'main iteration'
               REAL(KIND=dp) :: r_dec, r_inc   ! 'r_dec'=the decrease parameter and 'r_inc'=the increase parameter  
               REAL(KIND=dp) :: c              ! 'c'=the decrease parameter
               
               REAL(KIND=dp) :: step_tol       ! Step-length tolerance used in 'Clarke stationary' algorithm (same as the proximity measure)     
               REAL(KIND=dp) :: m_clarke       ! the descent parameter used in 'Clarke stationary' algorithm
           
               REAL(KIND=dp) :: f_0            ! the value of the objective function f=f_1-f_2 at the starting point x_0
               
               REAL(KIND=dp) :: f1_current, f2_current                                          ! the value of f_1 and f_2 at the current solution  
               REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: f1_current_block, f2_current_block   ! the value of f_1 and f_2 at the current solution  
               REAL(KIND=dp), DIMENSION(n_blocks) :: f1_0_block, f2_0_block                     ! the value of f_1 and f_2 at the current solution  
               REAL(KIND=dp), DIMENSION(n_blocks) :: f0_block                                   ! help variable  
               
               REAL(KIND=dp) :: f1_new, f2_new                                           ! the value of f_1 and f_2 at the new iteration point x_new (which is obtained from the previous main iteration)        
               REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: f1_new_block, f2_new_block    ! the value of f_1 and f_2 at the new iteration point x_new (which is obtained from the previous main iteration)        
               
               REAL(KIND=dp) :: change                           ! 'change' = f(x_new) - f(x_current)  (i.e. the change in the objective function value)
               REAL(KIND=dp), DIMENSION(n_blocks) :: change1     ! 'change1' = f1(x_new) - f1(x_current)  (i.e. the change in the DC component f1)
               REAL(KIND=dp), DIMENSION(n_blocks) :: change2     ! 'change2' = f2(x_new) - f2(x_current)  (i.e. the change in the DC component f2)
               
               REAL(KIND=dp), DIMENSION(nochange_mainit) :: fk_table  ! table for saving objective function values at the previous iteration points
               
               REAL(KIND=dp) :: norm                        ! the distance between two consecutive iteration points        
               REAL(KIND=dp) :: delta                       ! variable which can beused to relax the descent condition    
               
               REAL(KIND=dp) :: f_best                      ! help variable
    
               REAL(KIND=dp) :: start_time, middle_time, finish_time      ! start, middle and finish CPU time
               REAL(KIND=dp) :: start_time_main_it, finish_time_main_it   ! start and finish CPU time in one 'main iteration'
               
               REAL(KIND=dp) :: elapsed_time                  ! elapsed 'clock' time in seconds
               INTEGER :: clock_start, clock_end, clock_rate  ! start and finish 'clock' time   

               INTEGER :: size_b1 , size_b2    ! The biggest possible size of the bundles B_1 and B_2 
               
               INTEGER, DIMENSION(n_blocks) :: block_label   ! labels for each block
               INTEGER :: different_labels                   ! the number of different labels
               
               INTEGER :: bstart, bstart0  ! the start index of the block
               INTEGER :: bend, bend0      ! the end index of the block
               
               INTEGER  :: iter_counter        ! the number of main iterations executed 
               INTEGER  :: subprob_counter     ! the number of subproblems solved during the execution of the bundle algorithm
               INTEGER  :: stop_cond_counter   ! the number of times 'Clarke stationary' algorithm is used during the algorithm
               INTEGER  :: f_counter           ! the number of function values evaluated for a DC component during 
                                               ! the execution of the bundle algorithm (same for the DC component f_1 and f_2) (without Clarke stationary algorithm)

               INTEGER :: clarke_f_counter     ! the number of function values evaluated for f in 'Clarke stationary' algorithm during the execution of the bundle algorithm
               INTEGER :: clarke_sub_counter   ! the number of subgradients evaluated for f in 'Clarke stationary' algorithm during the execution of the bundle algorithm
               
               INTEGER  :: subgrad1_counter    ! the number of subgradients calculated for f_1 during the bundle algorithm (without Clarke stationary algorithm)
               INTEGER  :: subgrad2_counter    ! the number of subgradients calculated for f_2 during the bundle algorithm (without Clarke stationary algorithm)                 
               
               INTEGER :: help_mit_iter_counter    ! the number of iteration rounds executed in one 'main iteration'
               INTEGER :: help_subprob_counter     ! the number of subproblems solved in one 'main iteration'
               INTEGER :: help_f_counter           ! the number of function values evaluated for a DC component in one 'main iteration' (same for f_1 and f_2)
               INTEGER :: help_subgrad1_counter    ! the number of subgradients calculated for f_1 in one 'main iteration' (without Clarke stationary algorithm)
               INTEGER :: help_subgrad2_counter    ! the number of subgradients calculated for f_2 in one 'main iteration' (without Clarke stationary algorithm)
               INTEGER :: help_stop_cond_counter   ! the number of times CLarke stationarity was tested during the 'main iteration'
               
               INTEGER :: help_clarke_f_counter     ! the number of function values evaluated for f in 'Clarke stationary' algorithm
               INTEGER :: help_clarke_sub_counter   ! the number of subgradients evaluated for f in 'Clarke stationary' algorithm
               
               
               INTEGER :: reason_for_stop       ! the reason for stop during the 'main iteration'
                                                ! 0 - a new iteration point found
                                                ! 1 - stopping condition satisfied (Clarke stationarity)
                                                ! 2 - approximate stopping condition satisfied (i.e. the step-length beta* < step_tol)  
                                                ! 3 - the biggest possible number of rounds executed in the 'main iteration'
                                                ! 5 - the biggest possible number of rounds executed in the 'Clarke stationary' algorithm
               
               INTEGER :: max_threads           ! the maximum number of threads that can be used in parallellization
               INTEGER :: threads               ! the number of threads used in parallellization
               INTEGER :: max_sub_prob          ! the maximum number of subproblems solved at the same time
               
               INTEGER :: i, j, ind_best        ! help variables 
                                                    
               LOGICAL :: stop_alg              ! .TRUE. if the proximal double bundle algorithm DBDC can be terminated                                                        
               
               CALL cpu_time(start_time)                ! Start CPU timing     
               CALL cpu_time(middle_time)               ! Middle CPU timing     

               CALL SYSTEM_CLOCK(COUNT_RATE=clock_rate) ! Find the rate
               CALL SYSTEM_CLOCK(COUNT=clock_start)     ! Start timing 'clock' time            
 
!--------------------------------------------------------------------------------------
! The starting point for the test problem is selected from the 'starting_points' table using 'user_n' and 'start_ind' 
                
                ! The starting point of the first block is multiplied for other blocks
                DO i = 1,number_of_blocks
                   DO j = 1,size_of_block(1)
                       x_start(size_of_block(1)*(i-1)+j) = starting_points(start_ind,j)
                   END DO
                END DO

!--------------------------------------------------------------------------------------
                
                IF (boost == 1) THEN ! no boosting
                    DO i = 1,n_blocks
                        block_label(i) = i
                    END DO
                    different_labels = n_blocks
                    problem1 = problem10
                    problem2 = problem20
                    x_0 = x_start
                    size_blocks = size_blocks0
                ELSE ! boosting is used
                    x_0 = 0.0_dp
                    block_label = 0
                    problem1 = 0
                    problem2 = 0
                    size_blocks = 0
                    block_label(1) = 1
                    different_labels = 1
                    problem1(1) = problem10(1)
                    problem2(1) = problem20(1)
                    size_blocks(1) = size_blocks0(1)
                    DO i = 2,n_blocks
                        DO j = 1,i-1
                            IF (size_blocks0(i) == size_blocks0(j) .AND. &
                            & problem10(i) == problem10(j) .AND. problem20(i) == problem20(j)) THEN
                                block_label(i) = block_label(j)
                            END IF
                        END DO
                        IF (block_label(i) == 0) THEN
                            block_label(i) = different_labels + 1
                            different_labels = different_labels + 1
                            problem1(different_labels) = problem10(i)
                            problem2(different_labels) = problem20(i)
                            size_blocks(different_labels) = size_blocks0(i)
                        END IF
                    END DO
                    bstart = 1
                    DO i = 1,n_blocks
                        bend = bstart+size_blocks0(i)-1
                        f1_0_block(i) = f1_block(x_start(bstart:bend), problem10(i), size_blocks0(i), i)    ! the value of the DC component f_1 at x_0
                        f2_0_block(i) = f2_block(x_start(bstart:bend), problem20(i), size_blocks0(i), i)    ! the value of the DC component f_2 at x_0
                        f0_block(i) = f1_0_block(i) - f2_0_block(i)
                        bstart = bstart + size_blocks0(i)
                    END DO
                    DO i = 1,different_labels
                        f_best = 3.40282347*10.**38
                        ind_best = 0
                        DO j = 1,n_blocks
                            IF (block_label(j) == i) THEN
                                IF (f0_block(j) < f_best) THEN
                                    f_best = f0_block(j)
                                    ind_best = j
                                END IF
                            END IF
                        END DO
                        bstart0 = 1
                        DO j = 1,ind_best-1
                            bstart0 = bstart0 + size_blocks0(j)
                        END DO
                        bend0 = bstart0 + size_blocks0(ind_best) - 1
                        bstart = 1
                        DO j = 1,i-1
                            bstart = bstart + size_blocks(j)
                        END DO
                        bend = bstart + size_blocks(i) - 1
                        x_0(bstart:bend) = x_start(bstart0:bend0)
                    END DO
                END IF
                
                ALLOCATE(f1_current_block(different_labels))
                ALLOCATE(f2_current_block(different_labels))
                ALLOCATE(f1_new_block(different_labels))
                ALLOCATE(f2_new_block(different_labels))
                ALLOCATE(grad1(different_labels))
                ALLOCATE(grad2(different_labels))
                ALLOCATE(dd(different_labels))
                ALLOCATE(B1_table(different_labels))
                ALLOCATE(B2_table(different_labels))
               
                ! Initialization of the subgradients in the block form
                DO i=1, different_labels
                    ALLOCATE(grad1(i)%vector(size_blocks(i)))
                    ALLOCATE(grad2(i)%vector(size_blocks(i)))   
                    ALLOCATE(dd(i)%vector(size_blocks(i)))   
                END DO          
 
               
           !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^     
           ! <>  <>  <>  <>  <>  <>  <>  MAIN ITERATION STARTS  <>  <>  <>  <>  <>  <>  <>  <>
           !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^      
           
               !*************** PARAMETER VALUES FROM USER ARE CHECKED **********************
               !
               !     IF VALUES ARE INCORRECT THEN DEFAULT VALUES ARE USED FOR PARAMETERS
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
               
               ! ***** descent parameter 'r' *****
               IF ( (user_m<=0.0_dp) .OR. (user_m>=1.0_dp) ) THEN
                   m = 0.2_dp
               ELSE       
                   m = user_m 
               END IF
              
               ! ***** stopping tolerance 'crit_tol' *****
               IF (user_crit_tol <= 0.0_dp) THEN
                IF (optimality_condition == 1) THEN ! criticality
                   IF (user_n < 150) THEN 
                      crit_tol = user_n*0.005_dp
                   ELSE IF ( (150 <= user_n) .AND. (user_n <=200) ) THEN 
                      crit_tol = user_n*0.015_dp
                   ELSE IF (user_n > 200) THEN 
                      crit_tol = user_n*0.05_dp
                   END IF
                ELSE ! stationarity
                   IF (user_n < 50) THEN
                       crit_tol = (10.0_dp)**(-5)
                   ELSE IF ((50 <= user_n) .AND. (user_n <=200) ) THEN
                        crit_tol = (10.0_dp)**(-5)
                   ELSE IF (user_n > 200 ) THEN 
                        crit_tol = (10.0_dp)**(-4)                 
                   END IF
                END IF
               ELSE 
                   crit_tol = user_crit_tol
               END IF
               
               ! ***** proximity parameter 'eps' (for PBSDCcrit) *****
               IF (optimality_condition == 1) THEN !criticality
                   IF (user_eps <= 0.0_dp) THEN
                       eps = 0.1_dp
                   ELSE
                       eps = user_eps
                   END IF
               END IF
              
               ! ***** enlargement parameter 'eps' (for PBSDCstat) *****
               IF (optimality_condition == 2) THEN !stationarity
                   IF (user_eps_1 <= 0.0_dp) THEN
                       eps = 5*(10.0_dp)**(-5)
                   ELSE
                       eps = user_eps_1
                   END IF
               END IF
              
               ! ***** decrease parameter 'c1' *****
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

               ! ***** increase parameter 'C' *****
               IF ( user_r_inc <= 1.0_dp ) THEN
                   r_inc = (10.0_dp)**(7)
               ELSE
                   r_inc = user_r_inc
               END IF

               ! ***** bundle B1 size *****
               IF ( user_size_b1 <= 0 ) THEN
                   size_b1 = 0
                   DO i = 1,different_labels
                        IF (size_b1 < size_blocks(i)) THEN
                            size_b1 = size_blocks(i)
                        END IF
                   END DO
                   size_b1 = MIN(size_b1+5,100)
               ELSE
                   size_b1 = user_size_b1
               END IF

               ! ***** bundle B2 size *****
               IF ( user_size_b2 <= 0 ) THEN
                   size_b2 = 3
               ELSE
                   size_b2 = user_size_b2
               END IF      

               ! ***** decrease parameter 'c2' *****               
               IF ((user_c<=0.0_dp) .OR. (user_c>1.0_dp)) THEN
                   c = 0.1_dp
               ELSE
                   c = user_c
               END IF              
               
               !--------------------------------------------------------------------------
               !                      CLARKE STATIONARY ALGORITHM
               !--------------------------------------------------------------------------
               ! ***** descent parameter 'm_clarke' (for PBSDCstat) ***** 
               IF ((user_m_clarke<=0.0_dp) .OR. (user_m_clarke>=1.0_dp)) THEN
                   m_clarke = 0.01_dp
               ELSE
                   m_clarke = user_m_clarke
               END IF
               
               ! ***** step-lenght tolerance 'step_tol' (for PBSDCstat) ***** 
               IF (optimality_condition == 2) THEN !stationarity
                   IF (user_eps <= 0.0_dp) THEN
                       IF (user_n < 50) THEN
                           step_tol = 0.000001_dp
                       ELSE
                           step_tol = 0.00001_dp               
                       END IF
                   ELSE
                       step_tol = user_eps
                   END IF
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
               
               ! ***** maximum number of rounds 'mrounds_clarke' in the Clarke stationary algorithm *****
               IF ( mrounds_clarke <= 0 ) THEN
                   mrounds_clarke = 5000
               END IF                  
               
               ! ***** print option 'iprint' *****
               IF ( (iprint < -4 ) .OR. (iprint > 4) ) THEN   !executed if print value is wrong
                   iprint = 1
               END IF              
              !----------------------------------------------------------------------------
              
           !_______________________________________________________________________________
           !************************ STEP 0: PARAMETER INITIALIZATION *********************                
           
               x_current = x_0   ! the current iteration point is the starting point x_0
               
               bstart = 1 
               
               DO i = 1,different_labels
                   bend = bstart+size_blocks(i)-1
                   f1_current_block(i) = f1_block(x_0(bstart:bend), problem1(i), size_blocks(i), i)    ! the value of the DC component f_1 at x_0
                   f2_current_block(i) = f2_block(x_0(bstart:bend), problem2(i), size_blocks(i), i)    ! the value of the DC component f_2 at x_0

                   grad1(i)%vector = subgradient_f1_block(x_0(bstart:bend), problem1(i), size_blocks(i), i)    ! the subgradient of f_1 at x_0
                   grad2(i)%vector = subgradient_f2_block(x_0(bstart:bend), problem2(i), size_blocks(i), i)    ! the subgradient of f_2 at x_0
               
                   ! The bundles B_1 and B_2 are initialized
                   CALL init_bundle_b1(B1_table(i), size_b1, size_blocks(i))    
                   CALL init_bundle_b2(B2_table(i), size_b2, size_blocks(i))
                   
                   ! The first bundle element is added into the bundle B_1 and B_2 (i.e. the one corresponding to the starting point)
                   CALL add_first_element_b1(B1_table(i), grad1(i)%vector)
                   CALL add_first_element_b2(B2_table(i), grad2(i)%vector)
                   
                   bstart = bstart + size_blocks(i)
               END DO
               
               f_counter = 1                     ! one function value evaluated for both f_1 and f_2
               subgrad1_counter = 1              ! one subgradient calculated for f_1
               subgrad2_counter = 1              ! one subgradient calculated for f_1
               
               f1_current = 0.0_dp
               f2_current = 0.0_dp
               DO i = 1,n_blocks
                    f1_current = f1_current + f1_current_block(block_label(i))
                    f2_current = f2_current + f2_current_block(block_label(i))
               END DO
               f_0 = f1_current - f2_current     ! the value of the objective function at the starting point x_0 
               f0 = f_0                          ! the value at the starting point returned to the user
               
               iter_counter = 0             ! the number of 'main iterations' executed so far is zero
               subprob_counter = 0          ! the number of 'subproblems' solved so far is also zero
               stop_cond_counter = 0        ! the number of times approximate stopping condition is tested is zero
               stop_alg = .FALSE.           ! we cannot stop the proximal bundle algorithm
               
               clarke_f_counter = 0         ! the initialization of function value counter for 'Clarke stationary' algorithm
               clarke_sub_counter = 0       ! the initialization of subgradient counter for 'Clarke stationary' algorithm
               
               delta = 0.0_dp              
               
! --- --- --- Needed in OpenMP when we use PARALLELLIZATION --- --- ---   
               !$max_threads = omp_get_max_threads()
               !$max_sub_prob = give_max_size_b2(B2)+1
               !$threads = MIN(max_threads, max_sub_prob)
               !$CALL omp_set_num_threads(threads)    
! ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

              
               IF ( (ABS(iprint) >= 3) .AND. (ABS(iprint) <= 4) ) THEN        ! the basic/extended print of indermediate results     

                    WRITE(*,*) 'Main Iter:', 0, 'f(x):', f1_current - f2_current  

               END IF              
           !_______________________________________________________________________________
           !************************ STEP 0: END ******************************************    
           
           !----------------------------------------------------------------------------------   
           ! <>  <>  <>  <>  <>  <>  <> REPEATABLE PART BEGINS <>  <>  <>  <>  <>  <>  <>  <>
           !----------------------------------------------------------------------------------  
           
            DO WHILE ( ( .NOT. stop_alg) .AND. (iter_counter < mit ) .AND. &   ! This loop is repeated until the bundle algorithm DBDC can be terminated 
                        & ( (middle_time-start_time) < time_limit ) )          ! This test whether CPU time limit is met
            
                                                                         
               !_______________________________________________________________________________
               !************************ STEP 1: MAIN ITERATION ******************************* 

                iter_counter = iter_counter + 1     ! a new 'main iteration' is executed
                
                CALL cpu_time(start_time_main_it)
               
                ! The execution of main iteration
                CALL PBSDC_main_iteration( x_current, f1_current, f2_current, delta, user_n, &   
                            & x_new, f1_new, f2_new, f_0, B1_table, B2_table, iprint,  &
                            & crit_tol, eps, m, c, r_dec, r_inc, m_clarke, step_tol, mrounds,& 
                            & mrounds_clarke, reason_for_stop, help_mit_iter_counter, &
                            & help_subprob_counter, help_f_counter, help_subgrad1_counter, &
                            & help_subgrad2_counter, help_stop_cond_counter, &
                            & help_clarke_f_counter, help_clarke_sub_counter, &
                            & agg_used, stepsize_used, problem1, problem2, problem10, problem20, &
                            & n_blocks, size_blocks, size_blocks0, f1_new_block, f2_new_block, &
                            & f1_current_block, f2_current_block, L1, L2, optimality_condition, &
                            & block_label, different_labels, boost )  


                CALL cpu_time(finish_time_main_it)
                
                ! Different counters are updated
                subprob_counter = subprob_counter + help_subprob_counter
                f_counter = f_counter + help_f_counter
                stop_cond_counter = stop_cond_counter + help_stop_cond_counter
                subgrad1_counter = subgrad1_counter + help_subgrad1_counter
                subgrad2_counter = subgrad2_counter + help_subgrad2_counter 
                
                clarke_f_counter = clarke_f_counter + help_clarke_f_counter
                clarke_sub_counter = clarke_sub_counter + help_clarke_sub_counter

               !_______________________________________________________________________________
               !************************ STEP 1: END ****************************************** 
               
               IF (reason_for_stop == 0) THEN   ! in this case, a new iteration point is found in the 'main iteration'
               !_______________________________________________________________________________
               !************************ STEP 2: BUNDLE UPDATE*******************************
               
                    bstart = 1
                    DO i = 1,different_labels
                        bend = bstart+size_blocks(i)-1

                        ! Subgradients of f_1 and f_2 at the new point x_new
                        grad1(i)%vector = subgradient_f1_block(x_new(bstart:bend), problem1(i), size_blocks(i), i)
                        grad2(i)%vector = subgradient_f2_block(x_new(bstart:bend), problem2(i), size_blocks(i), i)
                        
                        change1(i) = f1_new_block(i) - f1_current_block(i)
                        change2(i) = f2_new_block(i) - f2_current_block(i)
                        
                        dd(i)%vector = x_new(bstart:bend) - x_current(bstart:bend)    ! the search direction
                        
                        CALL update_b1(B1_table(i), grad1(i)%vector, dd(i)%vector, change1(i))   ! bundle update for B_1
                        CALL update_b2(B2_table(i), grad2(i)%vector, dd(i)%vector, change2(i))   ! bundle update for B_2
                        
                        bstart = bstart + size_blocks(i)
                    END DO

                    change = f1_new - f2_new - (f1_current -f2_current) ! the change in the objective function value
                    
                    subgrad1_counter = subgrad1_counter + 1
                    subgrad2_counter = subgrad2_counter + 1
                    
                    IF (different_labels < n_blocks) THEN
                        bstart = 1
                        DO i=1,n_blocks
                            bstart0 = 1
                            DO j=1,block_label(i)-1
                                bstart0 = bstart0 + size_blocks(j)
                            END DO
                            bend0 = bstart0 + size_blocks(block_label(i)) - 1
                            bend = bstart+size_blocks0(i)-1
                            vect(bstart:bend) = x_new(bstart0:bend0) - x_current(bstart0:bend0) 
                            bstart = bstart + size_blocks0(i)
                        END DO
                    ELSE
                        vect = x_new - x_current
                    END IF
                    
                    norm = SQRT(DOT_PRODUCT(vect,vect))   ! distance between two consecutive iteration points
                                        
                    x_current = x_new                                   ! update of the current iteration point
                    f1_current = f1_new                                 ! update of the function value f_1
                    f2_current = f2_new                                 ! update of the function value f_2
                    f1_current_block = f1_new_block                     ! update of the function value f_1
                    f2_current_block = f2_new_block                     ! update of the function value f_2
                    
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
 
                    CALL cpu_time(middle_time)                          ! Middle CPU timing     
            
                    IF ( (ABS(iprint) >= 3) .AND. (ABS(iprint) <= 4) ) THEN        ! the basic/extended print of indermediate results
                        
                        WRITE(*,*) 'Main Iter:', iter_counter, 'f(x):', f1_current-f2_current, '  change:', change, 'norm: ', norm
                        IF (iprint > 3) THEN    ! IF iprint = 4 then the iteration point is printes
                            WRITE(*,*) 'The new iteration point:'
                            DO i = 1, user_n
                                WRITE(*,*) 'x(', i ,')*=',  x_current(i)
                            END DO
                            WRITE(*,*) ' '
                        END IF 
                        
                        IF ( ABS(iprint) == 4 ) THEN
                            WRITE(*,*) '--------------------------------------------------------------------------- '
                            WRITE(*,*) ' * * * * * DETAILIED INFORMATION ABOUT THE MAIN ITERATION', iter_counter ,' * * * * * ' 
                            WRITE(*,*) 'the number of rounds needed:', help_mit_iter_counter
                            WRITE(*,*) 'the number of subproblems solved:', help_subprob_counter
                            WRITE(*,*) 'the number of function values calculated for f_1:', help_f_counter
                            WRITE(*,*) 'the number of function values calculated for f_2:', help_f_counter
                            WRITE(*,*) 'the number of subgradients calculated for f_1:', help_subgrad1_counter
                            WRITE(*,*) 'the number of subgradients calculated for f_2:', help_subgrad2_counter
                            WRITE(*,*) 'the number of times Clarke stationary algorithm was tested:',help_stop_cond_counter 
                            
                            IF (help_stop_cond_counter > 0) THEN
                                WRITE(*,*) '--------------------------------------------------------------------------- '   
                                WRITE(*,*) ' * * * * * * * * DETAILS ABOUT CLARKE STATIONARY ALGORITHM * * * * * * * *  '   
                                WRITE(*,*) '--------------------------------------------------------------------------- '   
                                WRITE(*,*) 'The number of function values calculated for f', help_clarke_f_counter  
                                WRITE(*,*) 'The number of subgradients calculated for f', help_clarke_sub_counter   
                                WRITE(*,*) '--------------------------------------------------------------------------- '   
                            END IF
                            
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
                    
                    IF(ABS(iprint) == 4) THEN ! the extended print of indermediate results
                        WRITE(*,*) ' '
                        WRITE(*,*) '--------------------------------------------------------------------------- '
                        WRITE(*,*) ' * * * * * DETAILIED INFORMATION ABOUT THE MAIN ITERATION', iter_counter ,' * * * * * '
                        WRITE(*,*) 'the number of rounds needed:', help_mit_iter_counter
                        WRITE(*,*) 'the number of subproblems solved:', help_subprob_counter
                        WRITE(*,*) 'the number of function values calculated for f_1:', help_f_counter
                        WRITE(*,*) 'the number of function values calculated for f_2:', help_f_counter
                        WRITE(*,*) 'the number of subgradients calculated for f_1:', help_subgrad1_counter
                        WRITE(*,*) 'the number of subgradients calculated for f_2:', help_subgrad2_counter
                        WRITE(*,*) 'the number of times Clarke stationary algorithm was tested:',help_stop_cond_counter 
                            
                        IF (help_stop_cond_counter > 0) THEN    
                            WRITE(*,*) '--------------------------------------------------------------------------- '   
                            WRITE(*,*) ' * * * * * * * * DETAILS ABOUT CLARKE STATIONARY ALGORITHM * * * * * * * *  '   
                            WRITE(*,*) '--------------------------------------------------------------------------- '   
                            WRITE(*,*) 'The number of function values calculated for f', help_clarke_f_counter  
                            WRITE(*,*) 'The number of subgradients calculated for f', help_clarke_sub_counter   
                            WRITE(*,*) '--------------------------------------------------------------------------- '   
                        END IF                          
                            
                        WRITE(*,*) '--------------------------------------------------------------------------- '       
                        WRITE(*,*) 'CPU time used:', finish_time_main_it-start_time_main_it                             
                        WRITE(*,*) '--------------------------------------------------------------------------- '
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
           

           termination = reason_for_stop    ! the cause of the termination
           
           IF ( iter_counter >= mit ) THEN  ! the maximum number of 'main iterations' is executed and this causes the termination
               termination = 4 
           END IF
         
           IF ( middle_time-start_time > time_limit ) THEN  ! the CPU time limit is exceeded
               IF (optimality_condition == 2) THEN  ! approximate stationarity condition
                   termination = 6 
               ELSE ! approximate criticality condition
                   termination = 5
               END IF
           END IF
           
           IF (different_labels < n_blocks) THEN
                bstart = 1
                DO i=1,n_blocks
                    bstart0 = 1
                    DO j=1,block_label(i)-1
                        bstart0 = bstart0 + size_blocks(j)
                    END DO
                    bend0 = bstart0 + size_blocks(block_label(i)) - 1
                    bend = bstart+size_blocks0(i)-1
                    x_solution(bstart:bend) = x_current(bstart0:bend0) 
                    bstart = bstart + size_blocks0(i)
                END DO
           ELSE
                x_solution = x_current               ! the solution to the minimization problem
           END IF
           
           f_solution = f1_current - f2_current ! the objective function value at the solution
           
           ! the values of different counters
           counter(1) = iter_counter 
           counter(2) = subprob_counter
           counter(3) = f_counter
           counter(4) = subgrad1_counter
           counter(5) = subgrad2_counter 
           counter(6) = stop_cond_counter
           counter(7) = clarke_f_counter
           counter(8) = clarke_sub_counter
            
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
               WRITE(*,*) '---------------------------------------------------------------------------'
               IF (optimality_condition == 2) THEN
                   WRITE(*,*) 'The number of times Clarke stationary algorithm (CSA) was used:', stop_cond_counter  
                   WRITE(*,*) 'The number of function values computed for f in CSA:', clarke_f_counter
                   WRITE(*,*) 'The number of subgradients calculated for f in CSA: ', clarke_sub_counter
                   WRITE(*,*) '---------------------------------------------------------------------------'  
               ELSE     
                   WRITE(*,*) 'The number of times approximate criticality was tested:', stop_cond_counter  
                   WRITE(*,*) '---------------------------------------------------------------------------'                
               END IF
               WRITE(*,*) 'CPU time used:', finish_time-start_time                          
               WRITE(*,*) '--------------------------------------------------------------------------- '
               WRITE(*,*) 'Elapsed time:', elapsed_time
               WRITE(*,*) '---------------------------------------------------------------------------'                
           END IF 

           DEALLOCATE(f1_current_block)
           DEALLOCATE(f2_current_block)
           DEALLOCATE(f1_new_block)
           DEALLOCATE(f2_new_block)
           
           DO i = 1, different_labels
               CALL deallocate_b1(B1_table(i))       
               CALL deallocate_b2(B2_table(i)) 
               DEALLOCATE(grad1(i)%vector)
               DEALLOCATE(grad2(i)%vector)
               DEALLOCATE(dd(i)%vector)
           END DO
           
            DEALLOCATE(grad1)
            DEALLOCATE(grad2)
            DEALLOCATE(dd)
            DEALLOCATE(B1_table)
            DEALLOCATE(B2_table)  
           
           END SUBROUTINE PBSDC_algorithm     
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
           
            SUBROUTINE PBSDC_main_iteration( x_k, f1_k, f2_k, f_val, user_n, &
                            & x_new , f1_new, f2_new, f_0, B1_table, B2_table, iprint, &
                            & crit_tol, eps, m, c, r_dec, r_inc, m_clarke, step_tol, mrounds, &
                            & mrounds_clarke, reason_for_stop, iter_counter, &
                            & subprob_counter, fi_counter, subgrad1_counter, subgrad2_counter,&
                            & stop_cond_counter, clarke_f_counter, clarke_sub_counter, &
                            & agg_in_use, stepsize_used, problem1, problem2, problem10, problem20, &
                            & n_blocks, size_blocks, size_blocks0, f1_new_block, f2_new_block, &
                            & f1_k_block, f2_k_block, L1, L2, optimality_condition, &
                            & block_label, different_labels, boost )  
                            
            !
            ! Executes the 'main iteration' algorithm. It is needed to find the search direction at the current iteration point. 
            !
            ! INPUT: * 'x_k'               : The current iteration point
            !        * 'user_n'            : The dimension of the problem
            !        * 'n_blocks'          : The number of blocks
            !        * 'size_blocks'       : The sizes of truncated list of the blocks
            !        * 'size_blocks0'      : The original sizes of the blocks
            !        * 'problem1'          : The f1 of the problem
            !        * 'problem10'         : The f1 of the problem
            !        * 'problem2'          : The f2 of the problem
            !        * 'problem20'         : The f2 of the problem
            !        * 'f1_k' and 'f2_k'   : The values of DC components f_1 and f_2 at the current iteration point
            !        * 'f1_k_block' and 'f2_k_block'   : The values of DC components f_1 and f_2 at the current iteration point
            !        * 'f_val'             : The value of objective function used in descent condition
            !        * 'f_0'               : The value of the objective function f at a starting point x_0
            !        * 'L1'                : Overestimated Lipschitz constant of the DC component f1
            !        * 'L2'                : Overestimated Lipschitz constant of the DC component f1
            !        * 'iprint'            : Specifies print
            !        * 'block_label'       : Labels for each block
            !        * 'different_labels'  : The number of different labels
            !        * 'boost'             : Is boosting used, 1=no boosting, 2=boosting is used
            !
            !        * 'crit_tol'          : The stopping tolerance 
            !        * 'eps'               : The enlargement parameter
            !        * 'm'                 : The descent parameter used in main iteration
            !        * 'c'                 : The decrease parameter
            !        * 'r_dec' and 'r_inc' : The decrease and increase parameters   
            !        * 'm_clarke'          : The descent parameter used in 'Clarke stationary' algorithm 
            !        * 'step_tol'          : The step-length tolerance used in 'Clarke stationary' algorithm (i.e. proximity measure)
            !
            !        * 'mrounds'           : The maximum number of possible rounds in the 'main iteration'
            !        * 'mrounds_clarke'    : The maximum number of possible rounds in the 'Clarke stationary' algorithm
            !
            !        * 'agg_in_use'        : If .TRUE. then aggregation is used in the algorithm
            !        * 'stepsize_used'     : If .TRUE. then simple stepsize determination is used in the algorithm after each main iteration
            
            !        * 'optimality_condition'     : The used optimality condition
            !
            !
            ! OUTPUT: * 'x_new'               : the new iteration point
            !         * 'f1_new' and 'f2_new' : the new value of DC components f_1 and f_2 at 'x_new'
            !         * 'f1_new_block' and 'f2_new_block' : the new value of DC components f_1 and f_2 at 'x_new'
            !         * 'reason_for_stop'     : Indicates the cause of the termination in the 'main iteration' 
            !
            !         * 'iter_counter'       : the number of rounds needed in 'main iteration'
            !         * 'subprob_counter'    : the number of subproblems solved in 'main iteration'
            !         * 'fi_counter'         : the number of function values evaluated for a DC component in 'main iteration' (same for f_1 and f_2) (wihtout the ones used in Clarke stationary algorithm)
            !         * 'subgrad1_counter'   : the number of subgradients calculated for f_1 in 'main iteration' (without the ones used in Clarke stationary algorithm)
            !         * 'subgrad2_counter'   : the number of subgradients calculated for f_2 in 'main iteration' (without the ones used in Clarke stationary algorithm)
            !         * 'stop_cond_counter'  : the number of times approximate stopping condition is tested during the 'main iteration' 
            !         * 'clarke_f_counter'   : the number of function values evaluated for f in 'Clarke stationary' algorithm
            !         * 'clarke_sub_counter' : the number of subgradients evaluated for f in 'Clarke stationary' algorithm
            !
            ! INOUT: * 'B1_table' and B2_table'    : the bundles of the DC components f_1 and f_2
            !
            ! NOTICE: The dimensions of vectors 'x_k' and 'x_new' has to be same. (i.e. the dimensio of 'x_k' and 'x_new' has to be
            !         'user_n' when SUBROUTINE PBSDC_main_iteration is used in SUBROUTINE bundle_method.
            !
            !***********************************************************************************
               IMPLICIT NONE
            !**************************** NEEDED FROM THE USER *********************************  
            
               TYPE(kimppu1), DIMENSION(different_labels), INTENT(INOUT) :: B1_table      ! the bundle B_1 for the DC component f_1
               TYPE(kimppu2), DIMENSION(different_labels), INTENT(INOUT) :: B2_table      ! the bundle B_2 for the DC component f_2
               REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: x_k      ! the current iteration point
               REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: x_new    ! the new iteration point if 'reason_for_stop'=0 in the 'main iteration'

               REAL(KIND=dp), INTENT(IN) :: f1_k, f2_k        ! the value of f_1 and f_2 at x_k
               REAL(KIND=dp), DIMENSION(different_labels), INTENT(IN) :: f1_k_block, f2_k_block      ! the value of f_1 and f_2 at x_k
               REAL(KIND=dp), INTENT(IN) :: f_val             ! the value used in descent condition
               REAL(KIND=dp), INTENT(OUT) :: f1_new, f2_new   ! the value of f_1 and f_2 at x_new if 'reason_for_stop=0' in the 'main iteration'
               REAL(KIND=dp), DIMENSION(different_labels), INTENT(OUT) :: f1_new_block, f2_new_block  ! the value of f_1 and f_2 at x_new if 'reason_for_stop=0' in the 'main iteration'
               REAL(KIND=dp), INTENT(IN) :: f_0               ! the value of f at the starting point x_0
               
               REAL(KIND=dp), INTENT(IN) :: crit_tol       ! crit_tol=stopping tolerance 
               REAL(KIND=dp), INTENT(IN) :: eps            ! eps=enlargemant parameter
               REAL(KIND=dp), INTENT(IN) :: m              ! m=descent parameter
               REAL(KIND=dp), INTENT(IN) :: c              ! c=decrease parameter 
               REAL(KIND=dp), INTENT(IN) :: r_dec, r_inc   ! r_dec=decrease parameter and r_inc=increase parameter
               REAL(KIND=dp), INTENT(IN) :: m_clarke       ! the descent parameter used in 'Clarke stationary' algorithm 
               REAL(KIND=dp), INTENT(IN) :: step_tol       ! the step-length tolerance used in 'Clarke stationary' algorithm (i.e. proximity measure)
               
               INTEGER, INTENT(IN) :: different_labels                   ! the number of different labels
               INTEGER, INTENT(IN) :: boost                              ! is boosting used, 1=no boosting, 2=boosting is used 
               INTEGER, DIMENSION(n_blocks), INTENT(IN) :: block_label   ! labels for each block
               
               INTEGER, INTENT(IN) :: optimality_condition                      ! the used optimality condition
               INTEGER, INTENT(IN) :: user_n                                    ! the dimension of the problem
               INTEGER, INTENT(IN) :: n_blocks                                  ! the number of blocks
               INTEGER, DIMENSION(n_blocks), INTENT(IN) :: size_blocks          ! the sizes of truncated list of the blocks
               INTEGER, DIMENSION(n_blocks), INTENT(IN) :: size_blocks0         ! the original sizes of the blocks
               INTEGER, DIMENSION(n_blocks), INTENT(IN) :: problem1, problem10  ! the f1 of the problem
               INTEGER, DIMENSION(n_blocks), INTENT(IN) :: problem2, problem20  ! the f2 of the problem
               
               REAL(KIND=dp), INTENT(IN) :: L1      ! overestimated Lipschitz constant of the DC component f1
               REAL(KIND=dp), INTENT(IN) :: L2      ! overestimated Lipschitz constant of the DC component f2
               
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
               
               INTEGER :: mrounds                          ! the maximum number of possible rounds in the 'main iteration'
                                                           ! If mrounds<=0, then DEFAULT value 500 is used (This is also done in the SUBROUTINE bundle_method()
               INTEGER :: mrounds_clarke                   ! the maximum number of possible rounds in the 'Clarke statinonary' algorithm

               INTEGER, INTENT(OUT) :: reason_for_stop    ! When stationarity condition is used
                                                             ! 0 - a new iteration point found
                                                             ! 1 - the stopping condition is satisfied (i.e. Clarke stationarity)
                                                             ! 2 - the approximate stopping condition is satisfied (i.e. the step-length beta* < eps)
                                                             ! 3 - the maximum number 'mrounds' of rounds executed in one main iteration
                                                             ! 4 - the maximum number of 'main iterations' is executed  
                                                             ! 5 - the maximum number 'mrounds_clarke' of rounds executed in one 'Clarke stationary' alqorithm
                                                             ! 6 - the time limit is met during the execution of the alqorithm
                                                          ! When criticality condition is used
                                                             ! 0 - a new iteration point found
                                                             ! 1 - the stopping condition is satisfied (i.e. criticality)
                                                             ! 2 - the approximate stopping condition is satisfied (i.e. eps-criticality)
                                                             ! 3 - the maximum number 'mrounds' of rounds is executed in one main iteration
                                                             ! 4 - the maximum number of 'main iterations' is executed 
                                                             ! 5 - the execution time exceeded the time limit
               
               INTEGER, INTENT(OUT) :: iter_counter        ! the number of rounds used in 'main iteration'
               INTEGER, INTENT(OUT) :: subprob_counter     ! the number of subproblems solved in 'main iteration'
               INTEGER, INTENT(OUT) :: fi_counter          ! the number of function values evaluated for a DC component in 'main iteration' (same for f_1 and f_2) (without the ones used in Clarke stationary algorithm)
               INTEGER, INTENT(OUT) :: subgrad1_counter    ! the number of subgradients calculated for f_1 in 'main iteration' (without the ones used in Clarke stationary algorithm)
               INTEGER, INTENT(OUT) :: subgrad2_counter    ! the number of subgradients calculated for f_2 in 'main iteration' (without the ones used in Clarke stationary algorithm)
               INTEGER, INTENT(OUT) :: stop_cond_counter   ! the number of times Clarke stationary algortihm is entered during the 'main iteration'  

               INTEGER, INTENT(OUT) :: clarke_f_counter    ! the number of function values evaluated for f in 'Clarke stationary' algorithm
               INTEGER, INTENT(OUT) :: clarke_sub_counter  ! the number of subgradients evaluated for f in 'Clarke stationary' algorithm
               
               LOGICAL :: agg_in_use              ! If .TRUE. then aggregation is used in the algorithm 
               LOGICAL :: stepsize_used           ! If .TRUE. then step-size determination is used in the algorithm 
               

           !***************************** LOCAL VARIABLES ************************************
               
               TYPE(kimppu1) :: B1                                        ! the bundle B_1 for the DC component f_1
               TYPE(kimppu2) :: B2                                        ! the bundle B_2 for the DC component f_2
               REAL(KIND=dp), DIMENSION(user_n) :: d_t                    ! the search direction
               TYPE(block), DIMENSION(:), ALLOCATABLE :: d_t_block        ! the search direction
               REAL(KIND=dp), DIMENSION(user_n) :: y                      ! the new auxilary point
               REAL(KIND=dp), DIMENSION(user_n) :: new_grad1, new_grad2   ! the subgradients of f_1 and f_2 at y
               REAL(KIND=dp), DIMENSION(size_blocks(1)) :: agg_grad       ! the new aggregated subgradient of f_1
               REAL(KIND=dp), DIMENSION(user_n) :: helpgrad               ! 'help' subgradient
               REAL(KIND=dp), DIMENSION(size_blocks(1)) :: apu_grad       ! 'help' subgradient
               REAL(KIND=dp), DIMENSION(user_n) :: vect                   ! 'help' vector               
               REAL(KIND=dp), DIMENSION(user_n) :: help_d                 ! 'help' vector               
               REAL(KIND=dp), DIMENSION(n_blocks) :: f1_0_block           ! 'help' vector               
               REAL(KIND=dp), DIMENSION(n_blocks) :: f2_0_block           ! 'help' vector               
               REAL(KIND=dp), DIMENSION(n_blocks) :: f0_block             ! 'help' vector               
               TYPE(block), DIMENSION(:), ALLOCATABLE :: vect_block       ! 'help' vector               
               
               REAL(KIND=dp) :: norm1              ! ||\bxi_1(x_k)||
               REAL(KIND=dp) :: max_norm           ! the value of the maximum subgradient norm ||\bxi_{2,max}||
               REAL(KIND=dp) :: t                  ! the proximity parameter t
               REAL(KIND=dp) :: div_t              ! 1.0_dp divided by the proximity parameter t
               REAL(KIND=dp) :: t_min, t_max       ! the bounds for proximity parameter t
               
               REAL(KIND=dp) :: d_norm             ! the norm of the direction vector d_t
               REAL(KIND=dp) :: delta1, delta2     ! the predicted changes of f_1 and f_2, respectively
               REAL(KIND=dp) :: f1_y, f2_y         ! the function values of f_1 and f_2 at y
               REAL(KIND=dp), DIMENSION(different_labels) :: f1_y_block, f2_y_block         ! the function values of f_1 and f_2 at y
               REAL(KIND=dp) :: f_k                ! the function value of f at x_k
               REAL(KIND=dp) :: real_decrease      ! the real decrease of the objective function f
               REAL(KIND=dp) :: app_decrease       ! the approximation decrease of the objective function f
               REAL(KIND=dp) :: lin_err1, lin_err2 ! the linearization errors of f_1 and f_2 (calculated using y)
               REAL(KIND=dp) :: agg_lin_err        ! the new aggregated linearization error of f_1
               REAL(KIND=dp) :: error              ! help linearization error
               REAL(KIND=dp) :: theta          
               REAL(KIND=dp) :: eps1, tolerance     
               
               REAL(KIND=dp) :: help               ! 'help' variable  
               REAL(KIND=dp) :: f_best             ! 'help' variable  
               INTEGER :: ind_best                 ! 'help' variable
               
               REAL(KIND=dp) :: obj                     ! the value of the objective function at the solution 'lambda' (i.e. '|| bxi ||')
               REAL(KIND=dp), DIMENSION(user_n) :: bxi  ! the difference of the combined subgradients bxi = bxi_1^* - bxi_2^*

               INTEGER :: f_counter                ! help counter for the number of function values
               INTEGER :: sub_counter              ! help counter for the number of subgradients 
               INTEGER :: clarke_iter_counter      ! help counter for the number  of iterations in Clarke stationary algorithm 
               
               INTEGER :: s_counter                ! the number of subproblems solved during the current iteration round of the 'main iteration' algorithm
               INTEGER :: i, j, ind                ! help variables
               INTEGER :: bstart, bend             ! help variables
               INTEGER :: bstart0, bend0           ! help variables
               INTEGER :: N1, N2                   ! help variables
               INTEGER :: N1_before, N2_before     ! help variables
                           
               LOGICAL :: t_changed                ! .TRUE. if t has been changed during the previous round (.TRUE. also when initialization is done)
               LOGICAL :: was_b1_full              ! .TRUE. if during the previous round something was added to B_1 and B_1 was full before this insertion
               LOGICAL :: agg_used                 ! .TRUE. if aggregation element is in use
               LOGICAL :: stop_main_it             ! .TRUE. if the current 'main iteration' can be stopped
               
               LOGICAL :: stop_linesearch          ! .TRUE. if the current 'line_search' can be stopped
               REAL(KIND=dp) :: stepsize           ! stepsize determined during 'line search' into descent direction
               REAL(KIND=dp) :: test_f1            ! the value of f1 at the point tested in 'line search'
               REAL(KIND=dp) :: test_f2            ! the value of f2 at the point tested in 'line search'
               REAL(KIND=dp), DIMENSION(n_blocks) :: test_f1_block            ! the value of f2 at the point tested in 'line search'
               REAL(KIND=dp), DIMENSION(n_blocks) :: test_f2_block            ! the value of f2 at the point tested in 'line search'
               REAL(KIND=dp), DIMENSION(user_n) :: test_y  ! a point tested during 'line search'
                       

           !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^     
           ! <>  <>  <>  <>  <>  <>  <>  MAIN ITERATION STARTS  <>  <>  <>  <>  <>  <>  <>  <>
           !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
           
                ! Initializations
                ALLOCATE(vect_block(different_labels))
                ALLOCATE(d_t_block(different_labels))
                DO i=1,different_labels
                    ALLOCATE(vect_block(i)%vector(size_blocks(i)))  
                    ALLOCATE(d_t_block(i)%vector(size_blocks(i)))   
                END DO
                
                d_t = 0.0_dp
           
           !__________________________________________________________________________________     
           !************************ STEP 0: STOPPING CONDITION ******************************
                
               iter_counter = 1                   ! the first round is executed
               subprob_counter = 0                ! the number of subproblems solved
               fi_counter = 0                     ! the number of function values evaluated
               subgrad1_counter = 0               ! the number of subgradients calculated for f_1
               subgrad2_counter = 0               ! the number of subgradients calculated for f_2
               stop_cond_counter = 0              ! the number of times approximate stopping condition is tested
               clarke_f_counter = 0               ! the number of function values evaluated for f in 'Clarke stationary' algorithm
               clarke_sub_counter = 0             ! the number of subgradients evaluated for f in 'Clarke stationary' algorithm
                   
               agg_used = .FALSE.                 ! Aggregation is initialized to be .FALSE.
           
               stop_main_it = .FALSE.             ! We are NOT ready to stop
               
               bstart = 1
               DO i = 1,n_blocks
                   bend = bstart+size_blocks0(i)-1
                   vect(bstart:bend) = give_subgrad_b1(B1_table(block_label(i)),0) - give_subgrad_b2(B2_table(block_label(i)),0) ! vector: ' \bxi_1(x_k) - \bxi_2(x_k) '
                   bstart = bstart + size_blocks0(i)
               END DO
               
               IF (optimality_condition == 1) THEN !criticality
                   tolerance = 0.1_dp*crit_tol
               ELSE
                   tolerance = crit_tol
               END IF
               
               start: IF ( SQRT(DOT_PRODUCT(vect,vect)) < tolerance ) THEN  ! ||\bxi_1(x_k) - \bxi_2(x_k)|| <  crit_tol 
               !**>>**>>**>> APPROXIMATE CRITICALITY OBTAINED <<**<<**<<**
               
                   IF (optimality_condition == 1) THEN  ! Approximate criticality condition used
               
                       stop_main_it = .TRUE.  ! the 'main iteration' can be stopped
                       reason_for_stop = 1    ! the reason for stopping is criticality 
            
                       IF(iprint == 5) THEN  ! print everything        
                           WRITE(*,*) 'At STEP 0 the current point satisfies CRITICALITY'   
                       END IF
                       
                   ELSE  ! Approximate stationarity condition used
                   
                        d_t = 1.0_dp
                       
                        IF (different_labels < n_blocks) THEN
                            bstart = 1
                            DO i=1,n_blocks
                                bstart0 = 1
                                DO j=1,block_label(i)-1
                                    bstart0 = bstart0 + size_blocks(j)
                                END DO
                                bend0 = bstart0 + size_blocks(block_label(i)) - 1
                                bend = bstart+size_blocks0(i)-1
                                vect(bstart:bend) = x_k(bstart0:bend0) 
                                bstart = bstart + size_blocks0(i)
                            END DO
                        ELSE
                            vect = x_k    ! the solution to the minimization problem
                        END IF
                       
                       ! Clarke stationary algorithm is executed
                        CALL PBSDC_guaranteeing_clarke( vect, f1_k, f2_k, reason_for_stop, &
                                & x_new , f1_new, f2_new, d_t, &
                                & crit_tol, m_clarke, step_tol, mrounds_clarke, &
                                & clarke_iter_counter, f_counter, sub_counter, &
                                & problem1, problem2, user_n, &
                                & n_blocks, size_blocks, size_blocks0, f1_new_block, f2_new_block, &
                                & block_label, different_labels, boost )
                        
                        IF (different_labels < n_blocks) THEN
                            vect = x_new
                            x_new = 0.0_dp
                            bstart = 1
                            DO i = 1,n_blocks
                                bend = bstart+size_blocks0(i)-1
                                f1_0_block(i) = f1_block(vect(bstart:bend), problem10(i), size_blocks0(i), i)    ! the value of the DC component f_1 at x_0
                                f2_0_block(i) = f2_block(vect(bstart:bend), problem20(i), size_blocks0(i), i)    ! the value of the DC component f_2 at x_0
                                f0_block(i) = f1_0_block(i) - f2_0_block(i)
                                bstart = bstart + size_blocks0(i)
                            END DO
                            DO i = 1,different_labels
                                f_best = 10.0_dp**10
                                ind_best = 0
                                DO j = 1,n_blocks
                                    IF (block_label(j) == i) THEN
                                        IF (f0_block(j) < f_best) THEN
                                            f_best = f0_block(j)
                                            ind_best = j
                                        END IF
                                    END IF
                                END DO
                                bstart0 = 1
                                DO j = 1,ind_best-1
                                    bstart0 = bstart0 + size_blocks0(j)
                                END DO
                                bend0 = bstart0 + size_blocks0(ind_best) - 1
                                bstart = 1
                                DO j = 1,i-1
                                    bstart = bstart + size_blocks(j)
                                END DO
                                bend = bstart + size_blocks(i) - 1
                                x_new(bstart:bend) = vect(bstart0:bend0)
                            END DO
                        END IF
                        
                       ! Counters are updated   
                       clarke_f_counter = clarke_f_counter + f_counter                  
                       clarke_sub_counter = clarke_sub_counter + sub_counter                    
               
                       stop_cond_counter = stop_cond_counter + 1    ! the update of the stopping condition counter
                       stop_main_it = .TRUE.                        ! the 'main iteration' can be stopped

                   END IF
           !__________________________________________________________________________________         
           !************************ STEP 0: END *********************************************         
               ELSE start   
               !_______________________________________________________________________________
               !************************ STEP 1: PARAMETER INITIALIZATION *********************
           
                   f_k = f1_k - f2_k      ! the function value of f at x_k
                   
                   norm1 = 0.0_dp
                   max_norm = 0.0_dp
                   bstart = 1
                   DO i = 1,n_blocks
                       bend = bstart+size_blocks0(i)-1
                       vect(bstart:bend) = give_subgrad_b1(B1_table(block_label(i)),0) ! the vector \bxi_1(x_k)
                       bstart = bstart + size_blocks0(i)
                   END DO
                   norm1 = norm1 + SQRT(DOT_PRODUCT(vect,vect))  ! the value of the norm ||\bxi_1(x_k)||
                   
                   bstart = 1
                   DO i = 1,n_blocks
                       bend = bstart+size_blocks0(i)-1
                       vect(bstart:bend) = give_subgrad_b2(B2_table(block_label(i)),0) ! the vector \bxi_2(x_k)
                       bstart = bstart + size_blocks0(i)
                   END DO
                   max_norm = SQRT(DOT_PRODUCT(vect,vect))  ! the value of the norm ||\bxi_2(x_k)||
                   DO j = 1,give_size_b2(B2_table(1))
                       bstart = 1
                       DO i = 1,n_blocks
                           bend = bstart+size_blocks0(i)-1
                           vect(bstart:bend) = give_subgrad_b2(B2_table(block_label(i)),j) ! the vector \bxi_2
                           bstart = bstart + size_blocks0(i)
                       END DO
                       help = SQRT(DOT_PRODUCT(vect,vect))  ! the value of the norm ||\bxi_2||
                       IF (help > max_norm) THEN
                           max_norm = help
                       END IF
                   END DO

                   IF (optimality_condition == 1) THEN
                       eps1 = 0.5_dp * eps / ( MAX(L1,L2) )                   ! the value of epsilon_1
                       t_min = (0.5_dp * r_dec * eps1) / ( norm1 + max_norm ) ! the lower bound for t
                       t_max = r_inc * t_min                                  ! the upper bound for t
                       theta = r_dec * t_min * crit_tol                       ! theta
                   ELSE           
                       t_min = (0.5_dp * r_dec * eps) / ( norm1 + max_norm )  ! the lower bound for t
                       t_max = r_inc * t_min                                  ! the upper bound for t
                   END IF
                   
                   t = PBSDC_select_value_t(t_min,t_max)    ! the parameter t is selected         
            
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
                                    
                   IF (different_labels == 1 .AND. boost == 2) THEN                  
   
                       IF (agg_in_use) THEN                      ! 'agg_in_use' tells whether we use aggregation or not. If 'agg_in_use=.FALSE. then also 'agg_used'=.FALSE.
                           agg_used = is_agg_used(B1_table(1))   ! tells whether aggregation element is used or not? 
                       END IF   
                       
                   END IF
                   
                   IF ( agg_used ) THEN     ! IF we have .FALSE. here then the algorithm does NOT use AGGREGATION

                        bend = 1+size_blocks(1)-1
                        CALL PBSDC_subproblem_solver(x_k(1:bend), size_blocks(1), give_size_b1(B1_table(1))+2, &  ! subproblems are solved with aggregation
                                    &  B1_table(1), B2_table(1), t, s_counter)
                   
                   ELSE

                        bstart = 1
                        DO i = 1,different_labels
                            bend = bstart+size_blocks(i)-1
                            CALL PBSDC_subproblem_solver(x_k(bstart:bend), size_blocks(i), give_size_b1(B1_table(i))+1, &  ! subproblems are solved without aggregation
                                    &  B1_table(i), B2_table(i), t, s_counter)
                            bstart = bstart + size_blocks(i)          
                        END DO
                        
                   END IF   
                                                     
                   subprob_counter = subprob_counter + s_counter    ! the number of subproblems solved so far
                                   
                   bstart = 1
                   DO i = 1,different_labels
                       bend = bstart+size_blocks(i)-1
                       CALL add_glob_index(B2_table(i))                ! calculates the index of the subproblem yielding the global solution        
                       d_t(bstart:bend) = give_solution(B2_table(i))   ! the global solution d_t is selected  
                       d_t_block(i)%vector = d_t(bstart:bend)
                       bstart = bstart + size_blocks(i)                        
                   END DO
                   
                    IF (different_labels < n_blocks) THEN
                        bstart = 1
                        DO i=1,n_blocks
                            bstart0 = 1
                            DO j=1,block_label(i)-1
                                bstart0 = bstart0 + size_blocks(j)
                            END DO
                            bend0 = bstart0 + size_blocks(block_label(i)) - 1
                            bend = bstart+size_blocks0(i)-1
                            help_d(bstart:bend) = d_t(bstart0:bend0) 
                            bstart = bstart + size_blocks0(i)
                        END DO
                    ELSE
                        help_d = d_t        ! initial direction for stationarity check
                    END IF
                   
                   d_norm = SQRT(DOT_PRODUCT(help_d,help_d))  ! the norm ||d_t|| of the solution d_t is calculated
                   
                   ! Aggregated element is calculated. 
                   ! NOTICE: The aggregated element (if used) is now always added into bundle B_1 (also in those cases when the bundle B_1 is not full)

                   ! It is also possible to add the aggregated element into bundle B_1 only when the bundle B_1 is full
                   ! IF (is_full_b1(B1)) THEN

                   IF (different_labels == 1 .AND. boost == 2) THEN
                       IF (agg_in_use) THEN                     ! Is done only when 'agg_in_use'=.TRUE.
                           CALL add_glob_index(B2_table(1))
                           ind = give_solution_ind(B2_table(1))
                           apu_grad = give_subgrad_b2(B2_table(1),ind)
                           div_t = 1.0_dp / t
                           DO i = 1, give_n_b2(B2_table(1))
                              agg_grad(i) = (-div_t) * d_t_block(1)%vector(i) + apu_grad(i)
                           END DO
                           agg_lin_err = - give_decrease(B2_table(1))
                           agg_lin_err = agg_lin_err - (div_t) * DOT_PRODUCT(d_t_block(1)%vector,d_t_block(1)%vector)
                           agg_lin_err = agg_lin_err + give_linerr_b2(B2_table(1),ind)
                           CALL add_agg_element_b1(B1_table(1),agg_grad,agg_lin_err)
                       END IF
                   END IF

                   !END IF  
                                    
               !______________________________________________________________________________
               !************************ STEP 2: END *****************************************     
                   
                   IF (optimality_condition == 1) THEN
                      tolerance = theta
                   ELSE
                      tolerance = crit_tol
                   END IF
                   
               !->->->->->-> EXECUTION OF BRANCH BEGINS (2 POSSIBLE BRANCHES) <-<-<-<-<-<-<-<-
                   branches: IF (d_norm < tolerance) THEN                
               !->->->->->->->->->->->->->-> BRANCH 1 BEGIN <-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<- 

                      IF (optimality_condition == 1) THEN  ! Approximate criticality condition used
                      
                        !__________________________________________________________________________
                        !******************** STEP 3: APPROXIMATE STOPPING CONDITION **************
                      
                        N1_before = give_size_b1(B1_table(1)) + 1
                        N2_before = give_size_b2(B2_table(1)) + 1
                        
                        DO j = N1_before-1,1,-1
                            error = 0.0_dp
                            DO i = 1,n_blocks
                                error = error + give_linerr_b1(B1_table(block_label(i)), j)
                            END DO
                            IF (error > eps) THEN
                                ! delete the bundle element j
                                DO i = 1,different_labels
                                    CALL delete_b1_ind(B1_table(i), j)
                                END DO
                            END IF
                        END DO
                        
                        DO j = N2_before-1,1,-1
                            error = 0.0_dp
                            DO i = 1,n_blocks
                                error = error + give_linerr_b2(B2_table(block_label(i)), j)
                            END DO
                            IF (error > eps) THEN
                                ! delete the bundle element j
                                DO i = 1,different_labels
                                    CALL delete_b2_ind(B2_table(i), j)
                                END DO
                            END IF
                        END DO

                        N1 = give_size_b1(B1_table(1)) + 1  ! the current size of the bundle B_1 (NOTICE: the 'current bundle element' is also taken into account)
                        N2 = give_size_b2(B2_table(1)) + 1  ! the current size of the bundle B_2 (NOTICE: the 'current bundle element' is also taken into account)
                        
                        IF (different_labels > 1 .OR. boost == 1) THEN
                        
                            CALL init_bundle_b1(B1, N1, user_n)
                            CALL init_bundle_b2(B2, N2, user_n)
                            
                            DO j = 0, N1-1
                                error = 0.0_dp
                                bstart = 1
                                DO i = 1,n_blocks
                                    bend = bstart+size_blocks0(i)-1
                                    helpgrad(bstart:bend) = give_subgrad_b1(B1_table(block_label(i)), j)
                                    error = error + give_linerr_b1(B1_table(block_label(i)), j)
                                    bstart = bstart + size_blocks0(i)
                                END DO
                                IF (j == 0) THEN
                                    CALL add_first_element_b1(B1, helpgrad)
                                ELSE
                                    CALL add_element_b1(B1, helpgrad, error)
                                END IF
                            END DO
                                    
                            DO j = 0, N2-1
                                error = 0.0_dp
                                bstart = 1
                                DO i = 1,n_blocks
                                    bend = bstart+size_blocks0(i)-1
                                    helpgrad(bstart:bend) = give_subgrad_b2(B2_table(block_label(i)), j)
                                    error = error + give_linerr_b2(B2_table(block_label(i)), j)
                                    bstart = bstart + size_blocks0(i)
                                END DO
                                IF (j == 0) THEN
                                    CALL add_first_element_b2(B2, helpgrad)
                                ELSE
                                    CALL add_element_b2(B2, helpgrad, error)
                                END IF
                            END DO
                           
                            ! The norm minimization problem is solved. The solution is bxi = bxi_1^* - bxi_2^*
                            CALL norm_solver(bxi, obj, N1+N2, N1, N2, B1, B2, user_n)  

                            CALL deallocate_b1(B1)
                            CALL deallocate_b2(B2)
                            
                        ELSE ! all the blocks contain same number of variables (and same functions f1 and f2)

                            IF ( (N1_before == N1) .AND. (N2_before == N2)) THEN 
                                obj = d_norm / t
                            ELSE
                                ! The norm minimization problem is solved. The solution is bxi = bxi_1^* - bxi_2^*
                                CALL norm_solver(bxi(1:size_blocks(1)), obj, N1+N2, N1, N2, &
                                     & B1_table(1), B2_table(1), size_blocks(1))
                                
                                bstart = 1
                                DO i = 1,n_blocks
                                    bend = bstart+size_blocks0(i)-1
                                    helpgrad(bstart:bend) = bxi(1:size_blocks(1))
                                    bstart = bstart + size_blocks0(i)
                                END DO
                                
                                obj = SQRT(DOT_PRODUCT(helpgrad,helpgrad))
                            
                            END IF
                        
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
                           t = PBSDC_select_value_t(t_min,t_max)     ! the selection of the parameter t
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
                   
                      ELSE  ! Approximate stationarity condition used
                      !__________________________________________________________________________
                      !******************** STEP 3: CLARKE STATIONARITY CHECK *******************
                      
                        IF (different_labels < n_blocks) THEN
                            bstart = 1
                            DO i=1,n_blocks
                                bstart0 = 1
                                DO j=1,block_label(i)-1
                                    bstart0 = bstart0 + size_blocks(j)
                                END DO
                                bend0 = bstart0 + size_blocks(block_label(i)) - 1
                                bend = bstart+size_blocks0(i)-1
                                vect(bstart:bend) = x_k(bstart0:bend0) 
                                bstart = bstart + size_blocks0(i)
                            END DO
                        ELSE
                            vect = x_k     ! the point where stationarity is checked
                        END IF
                        
                        ! Clarke stationary algorithm is executed
                        CALL PBSDC_guaranteeing_clarke( vect, f1_k, f2_k, reason_for_stop, &
                            & x_new , f1_new, f2_new, help_d, &
                            & crit_tol, m_clarke, step_tol, mrounds_clarke, &
                            & clarke_iter_counter, f_counter, sub_counter, &
                            & problem1, problem2, user_n, &
                            & n_blocks, size_blocks, size_blocks0, f1_new_block, f2_new_block, &
                            & block_label, different_labels, boost )
                            
                        IF (different_labels < n_blocks) THEN
                            vect = x_new
                            x_new = 0.0_dp
                            bstart = 1
                            DO i = 1,n_blocks
                                bend = bstart+size_blocks0(i)-1
                                f1_0_block(i) = f1_block(vect(bstart:bend), problem10(i), size_blocks0(i), i)    ! the value of the DC component f_1 at x_0
                                f2_0_block(i) = f2_block(vect(bstart:bend), problem20(i), size_blocks0(i), i)    ! the value of the DC component f_2 at x_0
                                f0_block(i) = f1_0_block(i) - f2_0_block(i)
                                bstart = bstart + size_blocks0(i)
                            END DO
                            DO i = 1,different_labels
                                f_best = 10.0_dp**10
                                ind_best = 0
                                DO j = 1,n_blocks
                                    IF (block_label(j) == i) THEN
                                        IF (f0_block(j) < f_best) THEN
                                            f_best = f0_block(j)
                                            ind_best = j
                                        END IF
                                    END IF
                                END DO
                                bstart0 = 1
                                DO j = 1,ind_best-1
                                    bstart0 = bstart0 + size_blocks0(j)
                                END DO
                                bend0 = bstart0 + size_blocks0(ind_best) - 1
                                bstart = 1
                                DO j = 1,i-1
                                    bstart = bstart + size_blocks(j)
                                END DO
                                bend = bstart + size_blocks(i) - 1
                                x_new(bstart:bend) = vect(bstart0:bend0)
                            END DO
                        END IF
                       
                        !Counters are updated
                        clarke_f_counter = clarke_f_counter + f_counter                 
                        clarke_sub_counter = clarke_sub_counter + sub_counter   
                       
                        stop_cond_counter = stop_cond_counter + 1   ! the update of the stopping condition counter
                        
                        stop_main_it = .TRUE.                       ! Main iteration can be stopped
                        
                      END IF

                   !__________________________________________________________________________
                   !******************** STEP 3: APPROXIMATE STOPPING CONDITION **************  
                
               !->->->->->->->->->->->->->-> BRANCH 1 END <-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-
                        
                   ELSE branches
               !->->->->->->->->->->->->->-> BRANCH 2 BEGIN <-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-
                   !__________________________________________________________________________
                   !******************** STEP 4: DESCENT TEST ********************************

                       y = x_k + d_t      ! a new auxilary point y
                       
                       bstart = 1
                       f1_y = 0.0_dp
                       f2_y = 0.0_dp
                       DO i = 1,different_labels
                         bend = bstart+size_blocks(i)-1
                         f1_y_block(i) = f1_block(y(bstart:bend), problem1(i), size_blocks(i), i)      ! the value of f_1 at y
                         f2_y_block(i) = f2_block(y(bstart:bend), problem2(i), size_blocks(i), i)      ! the value of f_2 at y 
                         bstart = bstart + size_blocks(i)
                       END DO
                       DO i = 1,n_blocks
                         f1_y = f1_y + f1_y_block(block_label(i))
                         f2_y = f2_y + f2_y_block(block_label(i))
                       END DO

                       fi_counter = fi_counter + 1 ! one more objective function value calculated for a DC component
               
                       real_decrease = ( f1_y - f2_y ) -  f_k  ! the real decrease in the value of f    

                       app_decrease = 0.0_dp
                       DO i = 1,n_blocks
                          app_decrease = app_decrease + give_decrease(B2_table(block_label(i)))
                       END DO
                   
                       branch2: IF ( real_decrease <= (m * app_decrease + f_val ) ) THEN   ! the descent condition is dependent on the value 'f_val' 
                       !**>>**>>**>> NEW ITERATION POINT FOUND <<**<<**<<**  
                       
                           IF(stepsize_used) THEN 
                           !-**--**-  STEPSIZE DETERMINATION  -**--**-
                           IF (real_decrease < -0.1_dp) THEN
                              stepsize = 1.0_dp
                              stop_linesearch = .FALSE.
                              DO WHILE((.NOT. stop_linesearch))
                                 test_y = y + d_t
                                 test_f1 = 0.0_dp
                                 test_f2 = 0.0_dp
                                 bstart = 1
                                 DO i = 1,different_labels
                                    bend = bstart+size_blocks(i)-1
                                    test_f1_block(i) = f1_block(test_y(bstart:bend), problem1(i), size_blocks(i), i)
                                    test_f2_block(i) = f2_block(test_y(bstart:bend), problem2(i), size_blocks(i), i)
                                    bstart = bstart + size_blocks(i)
                                 END DO
                                 DO i = 1,n_blocks
                                    test_f1 = test_f1 + test_f1_block(block_label(i))     ! the value of f_1 at test_y
                                    test_f2 = test_f2 + test_f2_block(block_label(i))     ! the value of f_1 at test_y
                                 END DO
                                 fi_counter = fi_counter + 1                                 
                                 IF ( (test_f1 - test_f2 - f_k) <= ( 0.7_dp * m * app_decrease)) THEN 
                                    stepsize = stepsize + 1.0_dp
                                    y = test_y
                                    f1_y = test_f1                                 
                                    f2_y = test_f2  
                                    f1_y_block = test_f1_block
                                    f2_y_block = test_f2_block                                  
                                 ELSE
                                    stop_linesearch = .TRUE.  
                                 END IF
                              END DO
                           END IF 
                           !-**--**- STEPSIZE DETERMINATION END -**--**-
                           END IF
                           
                           x_new = y               ! the new iteration point is the current auxilary point
                           f1_new = f1_y           ! the value of f_1 at x_new
                           f2_new = f2_y           ! the value of f_2 at x_new
                           f1_new_block = f1_y_block           ! the value of f_1 at x_new
                           f2_new_block = f2_y_block           ! the value of f_2 at x_new
                           stop_main_it = .TRUE.   ! the 'main iteration' can be stopped
                           reason_for_stop = 0     ! the reason for stopping is the new iteration point
    
                           
                   !__________________________________________________________________________
                   !******************** STEP 4: END *****************************************         
                   
                       ELSE branch2 
                       !----------- SOME UPDATES IN VARIABLES BEGINS -------------------------  
                       
                   !__________________________________________________________________________     
                   !******************** STEP 5: BUNDLE UPDATE *******************************
                          IF (optimality_condition == 1) THEN
                             tolerance = eps1
                          ELSE
                             tolerance = eps
                          END IF
                          
                          update: IF ( ((f1_y - f2_y - f_0)>0 ) .AND. & 
                                                      & (d_norm > tolerance)) THEN 
                          !----------- PARAMETER t REDUCTION --------------------------------
                              t = t - r_dec * ( t - t_min )  ! the parameter t is updated
                              t_changed = .TRUE.             ! the parameter t was changed
 

                          ELSE update
                          !----------- BUNDLE AND PARAMETER UPDATE --------------------------
                                                  
                              delta2 = 0.0_dp
                              DO j = 1,different_labels  
                                i = give_solution_ind(B2_table(j))              ! the index of the subproblem yielding the global solution
                                vect_block(j)%vector = give_subgrad_b2(B2_table(j), i)   ! the subgradient of f_2 in the subproblem yielding the global solution
                                delta2 = delta2 - DOT_PRODUCT(vect_block(j)%vector, d_t_block(j)%vector) + &
                                         & give_linerr_b2(B2_table(j),i)  ! the value of delta_2
                              END DO
                            
                              delta1 = app_decrease - delta2                  ! the value of delta_1

                              IF (real_decrease > -m * app_decrease)  THEN    ! adjustment of the proximity parameter if necessary
                                  t =  t - c * ( t - t_min )
                              END IF
                              
                              bstart = 1
                              DO i = 1,different_labels
                                 bend = bstart+size_blocks(i)-1

                                 new_grad1(bstart:bend) = subgradient_f1_block(y(bstart:bend), problem1(i), size_blocks(i), i)   ! a subgradient of f_1 at y

                                 lin_err1 = f1_k_block(i) - f1_y_block(i) + &
                                            & DOT_PRODUCT(d_t_block(i)%vector, new_grad1(bstart:bend))  ! a new linearization error for f_1
                                                        
                                 !->>>>>>>>>>>> BUNDLE B_1 UPDATE BEGINS <<<<<<<<<<<<<<<<<<<<<<<<<<- 
                                 bundle1: IF (is_full_b1(B1_table(i))) THEN   
                                 ! bundle B_1 is full and a new element is added to the bundle B_1
                                 ! the overwritten element is written/taken down 
                                 
                                      CALL add_element_b1(B1_table(i), new_grad1(bstart:bend), lin_err1)
                                      was_b1_full = .TRUE.
                                
                                 ELSE bundle1
                                 ! bundle B_1 is NOT full and a new element is added to the bundle B_1

                                     CALL add_element_b1(B1_table(i), new_grad1(bstart:bend), lin_err1)
                                     was_b1_full = .FALSE.

                                 END IF bundle1  
                                 
                                 bstart = bstart + size_blocks(i)
                                 
                              END DO    
                              subgrad1_counter = subgrad1_counter + 1    ! a new subgradient was evaluated for f_1                            
                            !->>>>>>>>>>>> BUNDLE B_1 UPDATE ENDS <<<<<<<<<<<<<<<<<<<<<<<<<<<-


                            !->>>>>>>>>>>> BUNDLE B_2 AND PARAMETER UPDATE BEGINS <<<<<<<<<<<-                              

                            bundle2: IF( delta2 >= 0 .AND. (give_max_size_b2(B2_table(1))>0) ) THEN   ! bundle B_2 is updated (happens only when delta2 >= 0 and bundle size is larger than 1) (SOME OTHER INSERTION RULES ALSO POSSIBLE TO APPLY!)
                                
                                bstart = 1
                                DO i = 1,different_labels
                                  bend = bstart+size_blocks(i)-1
                                
                                  new_grad2(bstart:bend) = subgradient_f2_block(y(bstart:bend), problem2(i), size_blocks(i), i)  ! a subgradient of f_2 at y  
                                  lin_err2 = f2_k_block(i) - f2_y_block(i) + &
                                             & DOT_PRODUCT(d_t_block(i)%vector, new_grad2(bstart:bend))  ! a new linearization error for f_2  
                                
                                  CALL add_element_b2(B2_table(i), new_grad2(bstart:bend), lin_err2)   ! a new element is inserted into the bundle B_2
                                  !In the algorithm we never overwrite the 'bundle element' yielding the previous global solution of the search direction problem.
                                  
                                  bstart = bstart + size_blocks(i)
                                END DO  
                                subgrad2_counter = subgrad2_counter + 1              ! a new subgradient was evaluated for f_2
                   !__________________________________________________________________________     
                   !******************** STEP 5: END *****************************************
                   
                   !__________________________________________________________________________     
                   !******************** STEP 6: PARAMETER UPDATE ****************************

                                IF (different_labels < n_blocks) THEN
                                    bstart = 1
                                    DO i=1,n_blocks
                                        bstart0 = 1
                                        DO j=1,block_label(i)-1
                                            bstart0 = bstart0 + size_blocks(j)
                                        END DO
                                        bend0 = bstart0 + size_blocks(block_label(i)) - 1
                                        bend = bstart+size_blocks0(i)-1
                                        vect(bstart:bend) = new_grad2(bstart0:bend0) 
                                        bstart = bstart + size_blocks0(i)
                                    END DO
                                ELSE
                                    vect = new_grad2               ! the point where stationarity is checked
                                END IF
                                  
                                help = SQRT(DOT_PRODUCT(vect,vect)) ! the norm of the new subgradient of f_2 (subgradient is calculated 
                                                                                ! at the new auxilary point y)
                                IF(help > max_norm) THEN 
                                    max_norm = help            ! the updated value of the maximum subgradient norm in the bundle B_2
                                    IF (optimality_condition == 1) THEN
                                        t_min = (0.5_dp * r_dec * eps1) / ( norm1 + max_norm ) ! the updated value of the lower bound t_min
                                        theta = r_dec * t_min * crit_tol                       ! the updated value of the parameter theta
                                    ELSE
                                        t_min = (0.5_dp * r_dec * eps) / ( norm1 + max_norm ) ! the updated value of the lower bound t_min
                                    END IF
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
                  iter_counter = iter_counter -1 
              END IF       
           
              ! Deallocations
              DO i=1, different_labels
                DEALLOCATE(vect_block(i)%vector)  
                DEALLOCATE(d_t_block(i)%vector)  
              END DO
              DEALLOCATE(vect_block)
              DEALLOCATE(d_t_block)
              
           END SUBROUTINE PBSDC_main_iteration
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
           
           FUNCTION PBSDC_select_value_t(t_min,t_max) RESULT(t)
                IMPLICIT NONE
                REAL(KIND=dp), INTENT(IN) :: t_min, t_max
                REAL(KIND=dp) :: t

                t = t_min + 0.8_dp * (t_max-t_min) ! selection of t
                
                !OTHER OPTIONS:
                
!               t = 0.5_dp * (t_min + t_max)    ! selects the middle point from interval [t_min, t_max]
!               t = t_min + 0.8*(t_max - t_min)  

!               t = t_min                   ! selects the lower bound
!               t = t_max                   ! selects the upper bound
                
           END FUNCTION PBSDC_select_value_t
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
        !  |                      Guaranteeing Clarke stationarity                          |  |
        !  |                                                                                |  |
        !  ----------------------------------------------------------------------------------  |
        !***************************************************************************************        
           
            SUBROUTINE PBSDC_guaranteeing_clarke( x_k, f1_k, f2_k, reason_for_stop,&
                            & x_new , f1_new, f2_new, d_cause,&
                            & crit_tol, m, step_tol, mrounds, &
                            & iter_counter, f_counter, subgrad_counter, &
                            & problem1, problem2, user_n, &
                            & n_blocks, size_blocks, size_blocks0, f1_new_block, f2_new_block, &
                            & block_label, different_labels, boost )
                                                
            !
            ! Executes the 'Clarke stationary' algorithm. It is needed to guarantee Clarke stationarity at the current iteration point. 
            ! If the current point is not Clarke stationary, then this algorithm generates a descent direction.
            !
            ! INPUT: * 'x_k'               : The current iteratioin point
            !        * 'f1_k' and 'f2_k'   : The values of DC components f_1 and f_2 at the current iteration point
            !
            !        * 'd_cause'           : The search direction causing the execution of the Clarke stationary algorithm
            !        * 'crit_tol'          : The stopping tolerance
            !        * 'm'                 : The descent parameter
            !        * 'step_tol'          : The step-length tolerance (i.e. proximity measure)
            !
            !        * 'user_n'            : The dimension of the problem
            !        * 'n_blocks'          : The number of blocks
            !        * 'size_blocks'       : The sizes of truncated list of the blocks
            !        * 'size_blocks0'      : The original sizes of the blocks
            !        * 'problem1'          : The f1 of the problem
            !        * 'problem2'          : The f2 of the problem
            !        * 'mrounds'           : The maximum number of possible rounds in the 'Clarke stationary' algorithm
            
            !        * 'block_label'       : Labels for each block
            !        * 'different_labels'  : The number of different labels
            !        * 'boost'             : Is boosting used, 1=no boosting, 2=boosting is used 
            !
            ! OUTPUT: * 'x_new'              : the new iteration point (obtained if 'reason_for_stop' = 0)
            !         * 'f1_new' and 'f2_new': the new values of DC components f_1 and f_2 (obtained if 'reason_for_stop' = 0)
            !         * 'f1_new_block' and 'f2_new_block': the new values of DC components f_1 and f_2 (obtained if 'reason_for_stop' = 0)
            !
            !         * 'reason_for_stop'    : Indicates the cause of the termination in the 'Clarke stationarity' algorithm  
            !
            !         * 'iter_counter'       : the number of rounds needed in 'Clarke stationary' algorithm
            !         * 'f_counter'          : the number of function values evaluated for f in 'Clarke stationary' algorithm (same for f_1 and f_2)
            !         * 'subgrad_counter'    : the number of subgradients calculated for f in 'Clarke stationary' algorithm (same for f_1 and f_2)
            !
            ! NOTICE: The dimensions of vectors 'x_k' and 'x_new' has to be same. (i.e. the dimensio of 'x_k' and 'x_new' has to be
            !         'user_n' when SUBROUTINE PBSDC_guaranteeing_clarke is used in SUBROUTINE PBSDC_main_iteration.
            !
            !***********************************************************************************
               IMPLICIT NONE
            !**************************** NEEDED FROM THE USER ********************************* 
               REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: x_k      ! the current iteration point
               REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: d_cause  ! the search direction
               REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: x_new    ! the new iteration point if 'reason_for_stop'=0 in the 'Clarke stationary' algorithm

               REAL(KIND=dp), INTENT(IN)  :: f1_k, f2_k      ! the value of f_1 and f_2 at x_k
               REAL(KIND=dp), INTENT(OUT) :: f1_new, f2_new  ! the value of f_1 and f_2 at x_new if 'reason_for_stop=0' in the 'Clarke stationary' algorithm
               REAL(KIND=dp), DIMENSION(different_labels), INTENT(OUT) :: f1_new_block, f2_new_block  ! the value of f_1 and f_2 at x_new if 'reason_for_stop=0' in the 'Clarke stationary' algorithm
               
               REAL(KIND=dp), INTENT(IN) :: crit_tol      ! 'crit_tol'=stopping tolerance 
               REAL(KIND=dp), INTENT(IN) :: m             ! 'm'=descent parameter
               REAL(KIND=dp), INTENT(IN) :: step_tol      ! 'step_tol'=step-length tolerance (i.e. proximity measure)
               
               INTEGER, INTENT(IN):: user_n                               ! the dimension of the problem
               INTEGER, INTENT(IN):: n_blocks                             ! the number of blocks
               INTEGER, DIMENSION(n_blocks), INTENT(IN):: size_blocks     ! the sizes of truncated list of the blocks
               INTEGER, DIMENSION(n_blocks), INTENT(IN):: size_blocks0    ! the original sizes of the blocks
               INTEGER, DIMENSION(n_blocks), INTENT(IN):: problem1        ! the f1 of the problem
               INTEGER, DIMENSION(n_blocks), INTENT(IN):: problem2        ! the f2 of the problem
               
               INTEGER, INTENT(IN) :: different_labels                   ! the number of different labels
               INTEGER, INTENT(IN) :: boost                              ! is boosting used, 1=no boosting, 2=boosting is used 
               INTEGER, DIMENSION(n_blocks), INTENT(IN) :: block_label   ! labels for each block
               
               INTEGER :: mrounds                         ! the maximum number of possible rounds in the 'Clarke statinonary' algorithm
                                                          
               INTEGER, INTENT(OUT) :: reason_for_stop    ! 0 - a new iteration point found
                                                          ! 1 - stopping condition satisfied (Clarke stationarity)
                                                          ! 2 - approximate stopping condition satisfied (i.e. the step-length beta* < step_tol)
                                                          ! 5 - the maximum number of rounds executed in 'Clarke stationary' algorithm
                                                          
               INTEGER, INTENT(OUT) :: iter_counter       ! the number of rounds used in 'Clarke stationary' algorithm
               INTEGER, INTENT(OUT) :: f_counter          ! the number of function values evaluated for f in 'Clarke stationary' algorithm (same for f_1 and f_2)
               INTEGER, INTENT(OUT) :: subgrad_counter    ! the number of subgradients calculated for f in 'Clarke stationary' algorithm 
               

           !***************************** LOCAL VARIABLES ************************************
           
               TYPE(kimppu1) :: B                                ! the bundle B of f containing subgradients from '\partial f(x)'
               
               REAL(KIND=dp), DIMENSION(user_n) :: d                      ! the search direction
               REAL(KIND=dp), DIMENSION(user_n) :: x_d                    ! an auxiliary point
               REAL(KIND=dp), DIMENSION(user_n) :: new_grad               ! the subgradient of f at x
               REAL(KIND=dp), DIMENSION(user_n) :: new_grad1, new_grad2   ! the subgradients of f_1 and f_2 at x_d
               REAL(KIND=dp), DIMENSION(user_n) :: u_k                    ! the vector yielding minimum norm for quadratic norm minimization problem
               REAL(KIND=dp), DIMENSION(user_n) :: y                      ! a new auxiliary point            
               REAL(KIND=dp), DIMENSION(user_n) :: test_y                 ! another auxiliary point            
               
               REAL(KIND=dp) :: f_k              ! the value of the objective function at the current iteration point   
               REAL(KIND=dp) :: obj_u            ! the optimal value of the quadratin norm minimization problem 
               REAL(KIND=dp) :: norm_u           ! the norm for the vector u_k  
               REAL(KIND=dp) :: eps              ! stepsize used in subgradient calculation 
               REAL(KIND=dp) :: div_eps          ! 1.0_dp / eps     
               REAL(KIND=dp) :: real_decrease    ! the real decrease in the objective function f
               REAL(KIND=dp) :: f1_y, f2_y       ! the function values of f_1 and f_2 at y
               REAL(KIND=dp), DIMENSION(different_labels) :: f1_y_block, f2_y_block       ! the function values of f_1 and f_2 at y
               REAL(KIND=dp) :: test_f1, test_f2 ! the function values of f_1 and f_2 at test_y
               REAL(KIND=dp), DIMENSION(different_labels) :: test_f1_block, test_f2_block ! the function values of f_1 and f_2 at test_y
               REAL(KIND=dp) :: direc_der        ! directional derivative of f 
               REAL(KIND=dp) :: stepsize         ! stepsize 
               REAL(KIND=dp) :: descent_app      ! approximation of descent in the objective function value
               
               REAL(KIND=dp), ALLOCATABLE :: normu_table(:)  ! table for norm_u values
               INTEGER :: nochange_it         ! 0=not used, integer value 'k' = algorithm is stopped if after 'k' iterations change in the norm_u value is less than 'nochange_tol'
               REAL(KIND=dp) :: nochange_tol  ! tolerance when the change in norm_u is tested
               
               INTEGER :: bstart  ! the start index of the block
               INTEGER :: bend    ! the end index of the block
                            
               INTEGER :: size_b          ! the biggest possible bundle size for B
               INTEGER :: N               ! the current bundle size for B with the current element
               INTEGER :: i, k            ! help variable
                            
               LOGICAL :: stop_alg        ! .TRUE. if 'Clarke stationary' algorithm can be stopped
               LOGICAL :: stop_step_det   ! .TRUE. if the stepsize determination can be stopped
               

           !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^     
           !   <>  <>  <>  <>  <>  <>  <>  ALGORITHM STARTS  <>  <>  <>  <>  <>  <>  <>  <>
           !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
           !_______________________________________________________________________________
           !************************ STEP 0: INITIALIZATION *******************************    

            nochange_tol = 0.00000001_dp
            nochange_it = 20 ! 0 = not used
            ALLOCATE(normu_table(nochange_it))
                      
            IF (different_labels > 1 .OR. boost == 1) THEN            
            ! the size of the bundle B
               IF (user_size < 2) THEN
                   !size_b = user_n + 3 
                   size_b = MIN(user_n*2,200)
               ELSE
                   size_b = user_size
               END IF
            ELSE
               IF (user_size < 2) THEN
                   size_b = MIN(size_blocks(1)*2,200)
               ELSE
                   size_b = user_size
               END IF
            END IF
               
               eps = (10.0_dp)**(-3)     ! Stepsize used when subgradient is calculated
               div_eps = 1.0_dp / eps
           
            ! Initialization of iteration counters
               iter_counter = 1                       ! The first iteration round is executed
               f_counter = 0                          ! The number of function values evaluated for f
               subgrad_counter = 0                    ! The numeber of subgradients evaluated for f
               
               f_k = f1_k - f2_k                      ! The value of f at the current iteration point
               
               IF (SQRT(DOT_PRODUCT(d_cause, d_cause)) <= (10.0_dp)**(-10)) THEN
                   d = 1.0_dp
                   d = -d / SQRT(DOT_PRODUCT(d, d))    ! The search direction used to determine subgradient for f  
               ELSE
                   d = -d_cause / SQRT(DOT_PRODUCT(d_cause, d_cause))    ! The search direction used to determine subgradient for f   
               END IF    
               
            ! The bundles B is initialized         
            IF (different_labels > 1 .OR. boost == 1) THEN          
               CALL init_bundle_b1(B, size_b, user_n) 
            ELSE
               CALL init_bundle_b1(B, size_b, size_blocks(1)) 
            END IF
            ! Nothing is stored into bundle B      
               N = 0                        
               
               stop_alg = .FALSE.          ! Algorithm cannot be stopped
               

           !_______________________________________________________________________________    
           !******************** STEP 0: END **********************************************      
           
           
           !----------------------------------------------------------------------------------   
           ! <>  <>  <>  <>  <>  <>  <> REPEATABLE PART BEGINS <>  <>  <>  <>  <>  <>  <>  <>
           !----------------------------------------------------------------------------------         
           DO WHILE ( (.NOT. stop_alg) .AND. (iter_counter <= mrounds))                  
           !_______________________________________________________________________________
           !************************ STEP 1: NEW SUBGRADIENT ******************************         
           
               x_d = x_k + (eps) * d     ! Auxialiary point 'x_d' is calculated
               
               IF (different_labels > 1 .OR. boost == 1) THEN
                   bstart = 1
                   DO i = 1,n_blocks
                      bend = bstart+size_blocks0(i)-1
                      k = block_label(i)
                      new_grad1(bstart:bend) = subgradient_f1_block(x_d(bstart:bend), problem1(k), size_blocks(k), i)  ! Subgradient of f_1 at 'x_d'
                      new_grad2(bstart:bend) = subgradient_f2_block(x_d(bstart:bend), problem2(k), size_blocks(k), i)  ! Subgradient of f_2 at 'x_d'
                      bstart = bstart + size_blocks0(i)
                   END DO
               ELSE
                   new_grad1 = 0.0_dp
                   new_grad2 = 0.0_dp
                   bstart = 1
                   bend = bstart+size_blocks0(1)-1
                   new_grad1(bstart:bend) = subgradient_f1_block(x_d(bstart:bend), problem1(1), size_blocks(1), 1)  ! Subgradient of f_1 at 'x_d'
                   new_grad2(bstart:bend) = subgradient_f2_block(x_d(bstart:bend), problem2(1), size_blocks(1), 1)  ! Subgradient of f_2 at 'x_d'
               END IF
               
               ! new subgradient for f at x_k 
               new_grad = new_grad1 - new_grad2
               subgrad_counter = subgrad_counter +1 
               
               ! New subgradient is added into the bundle 'B'
               IF (N == 0) THEN 
                    IF (different_labels > 1 .OR. boost == 1) THEN
                      CALL add_first_element_b1(B, new_grad)  
                    ELSE
                      CALL add_first_element_b1(B, new_grad(1:size_blocks(1))) 
                    END IF
                    N = N + 1
               ELSE
                    IF (different_labels > 1 .OR. boost == 1) THEN
                      CALL add_element_b1(B, new_grad, 0.0_dp) 
                    ELSE
                      CALL add_element_b1(B, new_grad(1:size_blocks(1)), 0.0_dp) 
                    END IF
               END IF 
               
           !_______________________________________________________________________________    
           !******************** STEP 1: END **********************************************            
           
           !_______________________________________________________________________________
           !************************ STEP 2: CLARKE STATIONARITY **************************            
           
               N = give_size_b1(B) + 1      ! the current bundle size
               
               ! Minimum norm element 'u_k' is determined in the bundle 'B' (i.e. aggregated subgradient)
               IF (different_labels > 1 .OR. boost == 1) THEN
                   CALL PBSDC_quadratic_solver(x_k, user_n, N, B, 1.0_dp, u_k, obj_u) 
               ELSE
                   CALL PBSDC_quadratic_solver(x_k(1:size_blocks(1)), size_blocks(1), N, B, &
                          & 1.0_dp, u_k(1:size_blocks(1)), obj_u) 
                   bstart = 1
                   DO i=1,n_blocks
                        bend = bstart + size_blocks(1) - 1
                        u_k(bstart:bend) = u_k(1:size_blocks(1))
                        bstart = bend + 1
                   END DO
               END IF
               ! The solution u_k is has already negative sign. Thus d = u_k/norm_u
    
               norm_u = SQRT(DOT_PRODUCT(u_k,u_k))
               
                ! Check if norm_u has not changed in 'nochange_it' iterations
                IF (nochange_it > 0) THEN
                    IF (iter_counter < nochange_it+1) THEN
                        normu_table(iter_counter)=norm_u
                    ELSE
                        IF (ABS(normu_table(1)-norm_u) < nochange_tol) THEN
                            reason_for_stop = 8 ! norm_u is not changing
                            stop_alg = .TRUE.   ! the algorithm can be STOPPED
                        ELSE
                            DO i=1,nochange_it-1
                                normu_table(i)=normu_table(i+1)
                            END DO
                            normu_table(nochange_it)=norm_u
                        END IF
                    END IF
                END IF
            
               IF (norm_u < crit_tol) THEN  
                  
                  reason_for_stop = 1     ! Approximate Clarke stationarity is achieved
                  stop_alg = .TRUE.       ! Algorihtm can be stopped
                  CALL deallocate_b1(B)   ! Deallocate bundle B
                  
           !_______________________________________________________________________________    
           !******************** STEP 2: END **********************************************     
               
               ELSE            
           !_______________________________________________________________________________
           !************************ STEP 3: SEARCH DIRECTION *****************************            
               
               d =  u_k / norm_u          ! New search direction
               y = x_k + eps * d          ! New auxialiary point    
               
               f1_y = 0.0_dp
               f2_y = 0.0_dp

               IF (different_labels > 1 .OR. boost == 1) THEN
                   bstart = 1
                   DO i = 1,n_blocks
                      bend = bstart+size_blocks0(i)-1
                      k = block_label(i)
                      f1_y = f1_y + f1_block(y(bstart:bend), problem1(k), size_blocks(k), i)  ! f_1 at 'y'
                      f2_y = f2_y + f2_block(y(bstart:bend), problem2(k), size_blocks(k), i)  ! f_2 at 'y'
                      bstart = bstart + size_blocks0(i)
                   END DO
               ELSE
                   f1_y = n_blocks * f1_block(y(1:size_blocks(1)), problem1(1), size_blocks(1), 1)  ! f_1 at 'y'
                   f2_y = n_blocks * f2_block(y(1:size_blocks(1)), problem2(1), size_blocks(1), 1)  ! f_2 at 'y'
               END IF

               f_counter = f_counter +1 
               
               real_decrease = f1_y - f2_y - f_k      ! Real decraese in the objective function f

               direc_der = div_eps * real_decrease    ! Directional derivative at 'x_k' into the direction 'd'       
               
               IF (direc_der <= (-m * norm_u)) THEN   ! If .TRUE. decrease in the objective f is sufficient
                     reason_for_stop = 10             ! we will execute step-length determination at Step 4 
                     stop_alg = .TRUE.                ! algorithm can be stopped  
               END IF
               
               iter_counter = iter_counter +1         ! one more iteration is executed
               
               END IF 
               
           !_______________________________________________________________________________    
           !******************** STEP 3: END ********************************************** 
            END DO 
           !----------------------------------------------------------------------------------   
           !  <>  <>  <>  <>  <>  <>  <> REPEATABLE PART ENDS <>  <>  <>  <>  <>  <>  <>  <>
           !----------------------------------------------------------------------------------           
            
           !_______________________________________________________________________________
           !************************ STEP 4: STEP-LENGTH **********************************            
               IF ( (reason_for_stop > 1 .AND. (reason_for_stop < 8 .OR. reason_for_stop == 10)) &
                 & .AND. stop_alg ) THEN
           
               !-**--**- STEPSIZE DETERMINATION START -**--**-         
                  stepsize = 1.0_dp                      ! initial stepsize
                  stop_step_det = .FALSE.                ! stepsize determination cannot be stopped
                  DO WHILE((.NOT. stop_step_det))
                      test_y = x_k + stepsize * d        ! auxiliary test point
                      test_f1 = 0.0_dp
                      test_f2 = 0.0_dp
                      IF (different_labels > 1 .OR. boost == 1) THEN
                          bstart = 1
                          DO i = 1,different_labels
                              bend = bstart+size_blocks(i)-1
                              test_f1_block(i) = f1_block(test_y(bstart:bend), problem1(i), size_blocks(i), i)     ! the value of f_1 at test_y
                              test_f2_block(i) = f2_block(test_y(bstart:bend), problem2(i), size_blocks(i), i)     ! the value of f_2 at test_y
                              bstart = bstart + size_blocks(i) 
                          END DO
                          DO i = 1,n_blocks
                              test_f1 = test_f1 + test_f1_block(block_label(i))
                              test_f2 = test_f2 + test_f2_block(block_label(i))
                          END DO
                      ELSE
                          test_f1_block = f1_block(test_y(1:size_blocks(1)), problem1(1), size_blocks(1), 1)     ! the value of f_1 at test_y
                          test_f2_block = f2_block(test_y(1:size_blocks(1)), problem2(1), size_blocks(1), 1)     ! the value of f_2 at test_y
                          test_f1 = n_blocks * test_f1_block(1)
                          test_f2 = n_blocks * test_f2_block(1)
                      END IF
                      f_counter = f_counter + 1          
                      descent_app = -(10.0_dp)**(-4)     ! sufficient descent 
                      IF ( (test_f1 - test_f2 - f_k) > MIN(-m * stepsize * norm_u,descent_app )) THEN 
                         stepsize = 0.5_dp * stepsize    ! stepsize is decreased
                         IF (stepsize < step_tol) THEN  
                               stop_step_det = .TRUE.    
                         END IF
                      ELSE
                         stop_step_det = .TRUE.      
                         y = test_y
                         f1_y = test_f1
                         f2_y = test_f2
                         f1_y_block = test_f1_block
                         f2_y_block = test_f2_block
                      END IF
                  END DO
               !-**--**- STEPSIZE DETERMINATION END -**--**-  
               
                  CALL deallocate_b1(B)   ! Deallocate bundle B
                  
                  IF (stepsize >= step_tol) THEN    ! the step-length is big enough
                  
                     x_new = y                   ! the new iteration point is the auxilary point y
                     f1_new = f1_y               ! the value of f_1 at x_new
                     f2_new = f2_y               ! the value of f_2 at x_new
                     f1_new_block = f1_y_block   ! the value of f_1 at x_new
                     f2_new_block = f2_y_block   ! the value of f_2 at x_new
                     reason_for_stop = 0         ! the reason for stop is the new iteration point  
                     
                  ELSE
                  
                     reason_for_stop = 2     ! the reason for stop is that the approximate stopping condition is satisfied (i.e. the step-length beta* < step_tol)
                     
                  END IF 
               
                  iter_counter = iter_counter -1    ! one extra iteration in the counter needs to be removed 
                  
               END IF
           !_______________________________________________________________________________    
           !******************** STEP 4: END **********************************************             
           
              IF ( (.NOT. stop_alg ) ) THEN 
                  reason_for_stop = 5            ! the maximum number of rounds have been executed 
                               
              END IF       

              DEALLOCATE(normu_table)             
           
           !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^     
           !   <>  <>  <>  <>  <>  <>  <>  ALGORITHM ENDS  <>  <>  <>  <>  <>  <>  <>  <>
           !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^         
          
           
           END SUBROUTINE PBSDC_guaranteeing_clarke
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
        !  |                          QUADRATIC NORM SOLVER                                 |  |
        !  |                                                                                |  |        
        !  ----------------------------------------------------------------------------------  |
        !***************************************************************************************           

            SUBROUTINE PBSDC_quadratic_solver(x, NF, NA, B, t, d, obj) 
            ! Solves the quadratic norm problem in the 'Clarke stationary' algorithm
            ! Calls PLQDF1 by Ladislav Luksan
            !
            ! INPUT: * 'x'               : The current iteratioin point
            !        * 'NF' and 'NA'     : The number of variables and the size of the bundle 'B', respectively    
            !        * 'B'               : The bundle 'B'
            !        * 't'               : The proximity parameter
            !
            ! OUTPUT: * 'd'              : the vector giving minimum norm
            !         * 'obj'            : the optimal value of the quadratic norm problem
            !
            ! NOTICE: The dimensions of vectors 'x' and 'd' has to be same. (i.e. the dimensio of 'x' and 'd' has to be
            !         'user_n' when SUBROUTINE PBSDC_quadratic_solver is used in SUBROUTINE PBSDC_guaranteeing_clarke.
            
               IMPLICIT NONE               
            !**************************** NEEDED FROM USER ***********************************             
 
               TYPE(kimppu1), INTENT(IN) :: B                      ! the bundle B
               REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: x       ! the current iteration point
               REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: d       ! the vector giving minimum norm
               
               REAL(KIND=dp), INTENT(IN) :: t                      ! the proximity parameter t
               REAL(KIND=dp), INTENT(OUT) :: obj                   ! the optimal value of the quadratic norm problem
               
               INTEGER, INTENT(IN)  :: NF    ! the number of variables
               INTEGER, INTENT(IN)  :: NA    ! the bundle size of B is 'give_size_b1(B) + 1' 
               
              
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

               REAL(KIND=dp) :: u      ! help variable
               
               INTEGER :: j, k, l      ! help variables
               
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
               REAL(KIND=dp), DIMENSION(NF+1) :: G             ! Output of PLQDF1: gradient of the Lagrangian function
               REAL(KIND=dp), DIMENSION(NF*NA) :: AG           ! matrix whose columns are bundle subgradients 
               REAL(KIND=dp), DIMENSION((NF+1)*(NF+2)/2) :: AR ! Output of PLQDF1: triangular decomposition of kernel of the orthogonal projection
               REAL(KIND=dp), DIMENSION(NF*NC) :: CG           ! NOT significant since NC=0. matrix whose columns are normals of the linear constraints            
               
               ! .. Some other varibles ..    
               REAL(KIND=dp), DIMENSION(NF*NA) :: grad_m_b           ! subgradient matrix of B


           !************************** SUBPROBLEM SOLVER STARTS *********************************   
           
           !****************************** INITIALIZATIONS **************************************
           
               IX = 0             ! types of bounds: 0 - unbounded variables
               IA = 2             ! types of deviations
               
               u = 1.0_dp / t
               H = u              ! diagonal matrix
               
               grad_m_b = grad_matrix(B)  ! subgradient matrix of B 
               
               ! NO linearization errors
               AF = 0.0_dp      
  
               ! matrix containing subgradients   
               DO j = 1, NA
                   k = (j-1)*NF
                   DO l = 1, NF
                      AG(k+l) = grad_m_b(k+l) 
                   END DO
               END DO   
                        
               !Calls PLQDF1 by Ladislav Luksan
               CALL PLQDF1(NF,NA,NC,X,IX,XL,XU,AF,AFD,IA,IAA, &
                           & AG,AR,AZ,CF,IC,CL,CU,CG,G,H,S,MFP,KBF, &
                           & KBC,IDECF,ETA0,ETA2,ETA9,EPS7,EPS9,XNORM, &                       
                           & UMAX,GMAX,N,ITERQ)

               ! the vector giving minimum norm               
               DO j = 1, NF
                   d(j) = S(j) 
               END DO
              
              ! the optimal value of the quadratic norm problem
               obj = 0.5_dp * DOT_PRODUCT(d, d)         
    
                   
               
           END SUBROUTINE PBSDC_quadratic_solver
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

            SUBROUTINE PBSDC_subproblem_solver(x, NF, NA, B1, B2, t, subprob_counter) 
            ! Solves the subproblems in the original search direction problem (used in 'main iteration')
            ! Calls PLQDF1 by Ladislav Luksan
            !
            ! INPUT: * 'x'               : The current iteratioin point
            !        * 'NF' and 'NA'     : The number of variables and the size of the bundle 'B1', respectively    
            !        * 'B1' and 'B2'     : The bundle 'B1' and 'B2'
            !        * 't'               : The proximity parameter
            !
            ! OUTPUT: * 'subprob_counter'  : the number of subproblems solved
            !
            ! NOTICE: The dimensions of vectors 'x' has to be ' user_n' when SUBROUTINE PBSDC_subproblem_solver is used in SUBROUTINE PBSDC_main_iteration.
              
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
               REAL(KIND=dp) :: u, obj, a    ! help variables
               
               INTEGER :: i, j, k, l         ! help variables
               
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
                    
                       
               !->->->->->->->->->->->->->-> EACH SUBPROBLEM SOLVED <-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-
                    subproblems1: DO i = 0, subprob_counter   ! each subproblem is looked through
                    

                        grad_b2  = give_subgrad_b2(B2, i)  ! subgradient of B_2 in the subproblem i
                        alpha_b2 = give_linerr_b2(B2, i)   ! linearization error of B_2 in the subproblem i

                        
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
                              & AG,AR,AZ,CF,IC,CL,CU,CG,G,H,S,MFP,KBF, &
                              & KBC,IDECF,ETA0,ETA2,ETA9,EPS7,EPS9,XNORM, &                        
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
               !->->->->->->->->->->->->->-> EACH SUBPROBLEM SOLVED END <-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-                      
               !$OPM END PARALLEL DO

               subprob_counter = subprob_counter + 1 
               
               
           END SUBROUTINE PBSDC_subproblem_solver
        !.......................................................................................
        ! <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>        
        ! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  
        !| | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |        
        !** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** 
        !|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|



      END MODULE pbsdc