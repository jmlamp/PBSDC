        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !| |                                                                                  | |
        !| |                                                                                  | |
        !| |                                    PBSDC                                         | | 
        !| |     THE PROXIMAL BUNDLE METHOD FOR BLOCK-SEPARABLE NONSMOOTH DC OPTIMIZATION     | | 
        !| |                                 (version 1)                                      | |
        !| |                                                                                  | |
        !| |                                                                                  | |
        !| |          by Kaisa Joki and Jenni Lampainen (last modified September 2026)        | |
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
        !| |     The software is free for academic teaching and research purposes but we      | |
        !| |     ask you to refer the reference given below, if you use it.                   | |
        !| |                                                                                  | |
        !| |                                                                                  | |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*     
        !|                                                                                      |
        !|                                                                                      |
        !|    The code includes PBDC [2] and DBDC [3] methods by Kaisa Joki                     |
        !|    and both methods can be used in their original form in the code.                  |
        !|    Both methods are licensed by the MIT License.                                     |
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
        !|                                                                                      |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        !|                                                                                      |
        !|                                                                                      |
        !|   Codes include:                                                                     |
        !|                                                                                      |
        !|   tpbsdc.f95         - Main program for the new PBSDC method (this file)             |
        !|   constants.f95      - Double precision (also some parameters)                       |
        !|   bundle1.f95        - Bundle of DC component f_1                                    |
        !|   bundle2.f95        - Bundle of DC component f_2                                    |
        !|   functions.f95      - User-specified DC components f_1 and f_2 together with        |
        !|                        subgradients of DC components. Contains also user-specified   |
        !|                        initial values for parameters                                 |
        !|   norm_min.f95       - Solver for the norm minimization problem                      |
        !|   fun.f95            - Defines objective funtion and gradient of the norm            |
        !|                        minimization problem                                          |
        !|   pbsdc.f95          - New PBSDC method [1]                                          |
        !|                                                                                      |
        !|   pbdc.f95           - PBDC method by Kaisa Joki [2]                                 |
        !|   dbdc.f95           - DBDC method by Kaisa Joki [3]                                 |
        !|                                                                                      |
        !|   plqdf1.f           - Quadratic solver by Ladislav Luksan                           |
        !|   pvmm.f             - Variable metric method by Ladislav Luksan                     |
        !|   mqsubs.f           - Basic modules for PVMM by Ladislav Luksan                     |
        !|   pqsubs.f           - Matrix modules for PVMM by Ladislav Luksan                    |
        !|                                                                                      |
        !|   Makefile           - Makefile                                                      |
        !|                                                                                      |
        !|                                                                                      |
        !|                                                                                      |
        !|   To USE the software MODIFY tpbsdc.f95 and functions.f95 as needed                  |
        !|                                                                                      |
        !|                                                                                      |
        !|   Reference to PBSDC:                                                                |
        !|                                                                                      |
        !|   [1] J. Lampainen, K. Joki, A. M. Bagirov, S. Taheri and M. M. Mäkelä:              |
        !|       "Proximal bundle method for block-separable nonsmooth DC optimization".        |
        !|       Under review, (2026).                                                          |
        !|                                                                                      |
        !|   Reference to PBDC:                                                                 |
        !|                                                                                      |
        !|   [2] Kaisa Joki, Adil M. Bagirov, Napsu Karmitsa and Marko M. Mäkelä:               |
        !|       "A proximal bundle method for nonsmooth DC optimization utilizing              |
        !|       nonconvex cutting planes". J. Glob. Optim. 68(3), 501-535, (2017).             | 
        !|                                                                                      |
        !|   Reference to DBDC:                                                                 |
        !|                                                                                      |
        !|   [3] Kaisa Joki, Adil M. Bagirov, Napsu Karmitsa, Marko M. Mäkelä and Sona Taheri:  |
        !|       "Double bundle method for finding Clarke stationary points in nonsmooth        |
        !|       DC programming". SIAM J. Optim. 28(2), 1892-1919, (2018).                      |
        !|                                                                                      |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*       


      PROGRAM tpbsdc
      
         USE constants, ONLY : dp   ! double precision (i.e. accuracy)
         USE functions              ! INFORMATION from the USER
         USE bundle1                ! The BUNDLE of the DC component f_1
         USE bundle2                ! The BUNDLE of the DC component f_2
         USE dbdc                   ! DBDC method
         USE pbdc                   ! PBDC method
         USE pbsdc                  ! The new PBSDC method
         
        IMPLICIT NONE   
 
        REAL(KIND=dp), DIMENSION(10,30000)  :: starting_points       ! the starting point table
        REAL(KIND=dp), DIMENSION(10,30000)  :: starting_points_exp   ! the starting point table when exponent function involved
        
        INTEGER                           :: start_ind             ! the number of the used starting point
        
        CHARACTER(LEN=80)                 :: infile_starts         ! the file containing starting points
 
        REAL(KIND=dp) :: f_solution         ! The objective function value at the solution 'x_solution'
        REAL(KIND=dp) :: f0                 ! The objective function value at the starting point
        REAL(KIND=dp) :: f_solution_ave     ! The average objective function value at the starting point
        
        REAL(KIND=dp) :: time_limit         ! The CPU time limit for the problem
        REAL(KIND=dp) :: time               ! The CPU time
        REAL(KIND=dp) :: time_ave           ! The average CPU time
        
        REAL(KIND=dp) :: low_bound1, low_bound2, up_bound1, up_bound2    ! Bounds
        REAL(KIND=dp) :: c
        
        INTEGER, DIMENSION(8) :: counter    ! Contains the values of different counteres: 
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
                    
        INTEGER, DIMENSION(8) :: counter_ave    ! Contains the values of different counteres (averaged)
    
        INTEGER :: iprint                   ! Variable that specifies print option (specified by USER): 
                                            !   iprint = 0 : print is suppressed
                                            !   iprint = 1 : basic print of final result 
                                            !   iprint = -1: basic print of final result (without the solution vector)
                                            !   iprint = 2 : extended print of final result 
                                            !   iprint = -2: extended print of final result (without the solution vector)
                                            !   iprint = 3 : basic print of intermediate results and extended print of final results
                                            !   iprint = -3: basic print of intermediate results and extended print of final results (without the solution vector)
                                            !   iprint = 4 : extended print of intermediate results and extended print of final results 
                                            !   iprint = -4: extended print of intermediate results and extended print of final results (without the solution vectors)
                                            !
                                            ! If 'iprint' <= -5 .OR. 'iprint' >= 5 then DEFAULT value 'iprint'=1 is used    
                                            
        
        INTEGER, ALLOCATABLE :: problem1(:) ! problem number for f1 
        INTEGER, ALLOCATABLE :: problem2(:) ! problem number for f2 
        
        INTEGER :: n_problems               ! the number of problems (help variable)
        
        INTEGER :: group                    ! the solved group of problems
        INTEGER :: solver_ind               ! the solver used: 1=PBSDC, 2=DBDC, 3=PBDC
        INTEGER :: user_n                   ! the dimension of the problem
        
        INTEGER :: mit                      ! The maximum number of 'main iterations' (specified by USER).
                                            ! If 'mit' <=0 then DEFAULT value 'mit'=1000 is used

        INTEGER :: mrounds                  ! The maximum number of rounds during one 'main iteration' (specified by USER).
                                            ! If 'mrounds' <=0 then DEFAULT value 'mrounds'=5000 is used

        INTEGER :: mrounds_clarke           ! The maximum number of rounds during one 'Clarke stationarity' algorithm (specified by USER).
                                            ! If 'mrounds_clarke' <=0 then DEFAULT value 'mrounds_clarke'=5000 is used
                                            
        INTEGER :: termination              ! Tells the reason for termination in PBSDC  
                                            !-------------------------------------------------------------------------------------------------     
                                              ! When stationarity condition is used (same for PBSDCstat and DBDC)
                                                 ! 1 - the stopping condition is satisfied (i.e. Clarke stationarity)
                                                 ! 2 - the approximate stopping condition is satisfied (i.e. the step-length beta* < eps)
                                                 ! 3 - the maximum number 'mrounds' of rounds executed in one main iteration
                                                 ! 4 - the maximum number of 'main iterations' is executed  
                                                 ! 5 - the maximum number 'mrounds_clarke' of rounds executed in one 'Clarke stationary' alqorithm
                                                 ! 6 - the time limit is met during the execution of the alqorithm
                                                 ! 7 - after 'nochange_mainit' main iterations change in the objective value is less than 'nochange_tol'
                                            !-------------------------------------------------------------------------------------------------      
                                              ! When criticality condition is used (same for PBSDCcrit and PBDC)
                                                 ! 1 - the stopping condition is satisfied (i.e. criticality)
                                                 ! 2 - the approximate stopping condition is satisfied (i.e. eps-criticality)
                                                 ! 3 - the maximum number 'mrounds' of rounds is executed in one main iteration
                                                 ! 4 - the maximum number of 'main iterations' is executed 
                                                 ! 5 - the execution time exceeded the time limit
                                                 ! 7 - after 'nochange_mainit' main iterations change in the objective value is less than 'nochange_tol'

        LOGICAL :: agg_used                 ! If .TRUE. then aggregation is used in DBDC (specified by USER).
        LOGICAL :: stepsize_used            ! If .TRUE. then simple stepsize determination is done in DBDC after each 'main iteration' (specified by USER).

        INTEGER :: dim_loop                      ! Defines the dimension(s) in the considered problem 
                                                 !    1 - the dimensions are 2, 5, 10, 50, 100, 200
                                                 !    2 - the dimensions are 10, 25, 50, 100, 250, 500
                                                 
        INTEGER :: optimality_condition     ! 1=approximate criticality condition, 2=approximate stationarity condition
                                                 
        INTEGER :: ind_start                                         
        INTEGER :: ind_finish   

        INTEGER :: nochange_mainit       ! 0=not used, integer value 'k' = algorithm is stopped if after 'k' main iterations change in the objective value is less than 'nochange_tol'
        REAL(KIND=dp) :: nochange_tol    ! tolerance for change in the objective value

        INTEGER :: block_size_label

        INTEGER :: ind_start_block_size                                      
        INTEGER :: ind_finish_block_size
        
        INTEGER :: block_number_label
        
        INTEGER :: ind_start_block_number                                    
        INTEGER :: ind_finish_block_number   

        INTEGER :: boost                 ! 1 = no boosting, 2 = boosting is used     

        INTEGER :: number_of_blocks0                ! the number of blocks in the block-separable function
        INTEGER, ALLOCATABLE :: size_of_block0(:)   ! the number of variables in each block    

        INTEGER :: bound_label1, bound_label2       ! bound labels
        
        INTEGER :: i, j, k, l, kmax
        
        LOGICAL :: start_normal
        
        CHARACTER(LEN=20) :: solver_name

        CHARACTER*80 outfi0
        CHARACTER*80 outfi1
        CHARACTER*80 outfi2
        CHARACTER*80 outfi3
        CHARACTER*80 outfi0_ave
        CHARACTER*80 outfi1_ave
        CHARACTER*80 outfi2_ave
        CHARACTER*80 outfi3_ave
         
      !--------------------------------
      ! ** start PBSDC parameters ** 
      !--------------------------------
        
        time_limit = 900           ! 15 min time limit
        
        ! Default for these three is 10000      
        mrounds = 10000           ! maximum number of rounds during one 'main iteration'
        mit = 100000              ! maximum number of 'main iteration'
        mrounds_clarke = 10000    ! maximum number of rounds during one 'Clarke stationary' algorithm
          
        iprint = 1                ! basic print of intermediate results and extended print of final results
          
        agg_used = .TRUE.         ! Aggregation is used (default=TRUE)
        stepsize_used = .FALSE.   ! Simple stepsize determination is not used  

        boost = 2                 ! 1 = no boosting, 2 = boosting is used (default)
        
        !group = 0      ! The solved group of problems: separable problems
        group = 1      ! The solved group of problems: block-separable problems
        
        solver_ind = 1    ! PBSDC
        !solver_ind = 2   ! DBDC
        !solver_ind = 3   ! PBDC
        
        optimality_condition = 1    ! approximate criticality condition
        !optimality_condition = 2    ! approximate stationarity condition
        
        nochange_mainit = 0        ! 0=not used, integer value 'k' = algorithm is stopped if after 'k' main iterations change in the objective value is less than 'nochange_tol'
        nochange_tol = 0.00001_dp  ! tolerance for change in the objective value
        
        dim_loop = 1              ! 1 = problems in group 0
        block_size_label = 1      ! 1 = small problems, 2 = large problems (group 1)
        block_number_label = 1    ! 1 = small problems, 2 = large problems (group 1)
        
      !-----------------------------
      ! ** end PBSDC parameters **
      !-----------------------------
      
      ! Starting points are read from the files
        
        infile_starts = 'startpoints10.txt'
        
        OPEN(78,file=infile_starts,status='old',form='formatted')
            DO k=1,10
               READ(78,*) (starting_points(k,j),j=1,30000)
            END DO
        CLOSE(78)    

        infile_starts = 'startpoints5.txt'
        
        OPEN(78,file=infile_starts,status='old',form='formatted')
            DO k=1,10
               READ(78,*) (starting_points_exp(k,j),j=1,30000)
            END DO
        CLOSE(78) 
        
        SELECT CASE(solver_ind)
            CASE(1)
               solver_name = 'PBSDC'
               IF (boost == 2) THEN ! boosting is used
                   outfi0 = 'results_PBSDC_boost_group0.txt'
                   outfi1 = 'results_PBSDC_boost_group1.txt'
                   outfi2 = 'results_PBSDC_boost_group2.txt'
                   outfi3 = 'results_PBSDC_boost_group3.txt'
                   outfi0_ave = 'results_PBSDC_boost_group0_ave.txt'
                   outfi1_ave = 'results_PBSDC_boost_group1_ave.txt'
                   outfi2_ave = 'results_PBSDC_boost_group2_ave.txt'
                   outfi3_ave = 'results_PBSDC_boost_group3_ave.txt'
               ELSE ! no boosting
                   outfi0 = 'results_PBSDC_group0.txt'
                   outfi1 = 'results_PBSDC_group1.txt'
                   outfi2 = 'results_PBSDC_group2.txt'
                   outfi3 = 'results_PBSDC_group3.txt'
                   outfi0_ave = 'results_PBSDC_group0_ave.txt'
                   outfi1_ave = 'results_PBSDC_group1_ave.txt'
                   outfi2_ave = 'results_PBSDC_group2_ave.txt'
                   outfi3_ave = 'results_PBSDC_group3_ave.txt'
               END IF
            CASE(2)
               solver_name = 'DBDC'
               outfi0 = 'results_DBDC_group0.txt'
               outfi1 = 'results_DBDC_group1.txt'
               outfi2 = 'results_DBDC_group2.txt'
               outfi3 = 'results_DBDC_group3.txt'
               outfi0_ave = 'results_DBDC_group0_ave.txt'
               outfi1_ave = 'results_DBDC_group1_ave.txt'
               outfi2_ave = 'results_DBDC_group2_ave.txt'
               outfi3_ave = 'results_DBDC_group3_ave.txt'
            CASE(3)
               solver_name = 'PBDC'
               outfi0 = 'results_PBDC_group0.txt'
               outfi1 = 'results_PBDC_group1.txt'
               outfi2 = 'results_PBDC_group2.txt'
               outfi3 = 'results_PBDC_group3.txt'
               outfi0_ave = 'results_PBDC_group0_ave.txt'
               outfi1_ave = 'results_PBDC_group1_ave.txt'
               outfi2_ave = 'results_PBDC_group2_ave.txt'
               outfi3_ave = 'results_PBDC_group3_ave.txt'
        END SELECT
        
        !-------------------------------------------------------------------         
        !                             PRINT
        !-------------------------------------------------------------------    
        SELECT CASE(group)   
          CASE(0)
            OPEN(45,file=outfi0) 
            WRITE(45,*)  'Results obtained for Group 0 with ', solver_name
            OPEN(55,file=outfi0_ave) 
            WRITE(55,*)  'Results obtained for Group 0 with ', solver_name
            
          CASE(1)   
            OPEN(45,file=outfi1)  
            WRITE(45,*)  'Results obtained for Group 1 with ', solver_name
            OPEN(55,file=outfi1_ave) 
            WRITE(55,*)  'Results obtained for Group 1 with ', solver_name
        END SELECT
        
        IF (solver_ind == 1) THEN
            WRITE(45,*)  'Option for boost ', boost, '1=no boosting, 2=boosting is used'
            WRITE(55,*)  'Option for boost ', boost, '1=no boosting, 2=boosting is used'
        END IF
        
        WRITE(45,*)  ' '            
        WRITE(45,*)  'Parameters are defaults'
        WRITE(45,*)  ' '      
        WRITE(45,*)  'Maximum number of rounds during one main iterations =', mrounds     
        WRITE(45,*)  'Maximum number of main iterations =', mit   
        WRITE(45,*)  ' '      
        WRITE(45,*)  'Number of main iterations used between when objective values are checked =', nochange_mainit     
        WRITE(45,*)  'Tolerance for the change in the objective value =', nochange_tol

        IF (solver_ind == 2 .OR. (solver_ind == 1 .AND. optimality_condition == 2)) THEN            
            WRITE(45,*)  'Maximum number of  rounds during one Clarke stationary algorithm =', mrounds_clarke     
        END IF
            
        WRITE(45,*)  ' ' 
        WRITE(45,*)  'Termination explanations ' 
            
        SELECT CASE(solver_ind)
        
            CASE(1) ! PBSDC
                IF (optimality_condition == 2) THEN  ! approximate stationarity condition
                    WRITE(45,*)  ' 1 - the stopping condition is satisfied (i.e. Clarke stationarity)'
                    WRITE(45,*)  ' 2 - the approximate stopping condition is satisfied (i.e. the step-length beta* < user_eps)'
                    WRITE(45,*)  ' 3 - the maximum number mrounds of rounds is executed in one main iteration'
                    WRITE(45,*)  ' 4 - the maximum number of main iterations is executed'
                    WRITE(45,*)  ' 5 - the maximum number mrounds_clarke of rounds is executed in one Clarke stationary algorithm'
                    WRITE(45,*)  ' 6 - the time limit is met during the execution of the alqorithm' 
                    WRITE(45,*)  ' 7 - after "nochange_mainit" main iterations change in the objective is less than "nochange_tol"' 
                    WRITE(45,*)  ' 8 - if norm_u does not change in Clarke stationarity algorithm' 
                ELSE  ! approximate criticality condition
                    WRITE(45,*)  ' 1 - the stopping condition is satisfied (i.e. criticality)'
                    WRITE(45,*)  ' 2 - the approximate stopping condition is satisfied (i.e. eps-criticality)'
                    WRITE(45,*)  ' 3 - the maximum number mrounds of rounds is executed in one main iteration'
                    WRITE(45,*)  ' 4 - the maximum number of main iterations is executed'
                    WRITE(45,*)  ' 5 - the execution time exceeded the time limit'
                    WRITE(45,*)  ' 7 - after "nochange_mainit" main iterations change in the objective is less than "nochange_tol"' 
                END IF
            
            CASE(2) ! DBDC        
                WRITE(45,*)  ' 1 - the stopping condition is satisfied (i.e. Clarke stationarity)'
                WRITE(45,*)  ' 2 - the approximate stopping condition is satisfied (i.e. the step-length beta* < user_eps)'
                WRITE(45,*)  ' 3 - the maximum number mrounds of rounds is executed in one main iteration'
                WRITE(45,*)  ' 4 - the maximum number of main iterations is executed'
                WRITE(45,*)  ' 5 - the maximum number mrounds_clarke of rounds is executed in one Clarke stationary algorithm'
                WRITE(45,*)  ' 6 - the time limit is met during the execution of the alqorithm'   
                WRITE(45,*)  ' 7 - after "nochange_mainit" main iterations change in the objective is less than "nochange_tol"' 
            
            CASE(3) ! PBDC
                WRITE(45,*)  ' 1 - the stopping condition is satisfied (i.e. criticality)'
                WRITE(45,*)  ' 2 - the approximate stopping condition is satisfied (i.e. eps-criticality)'
                WRITE(45,*)  ' 3 - the maximum number mrounds of rounds is executed in one main iteration'
                WRITE(45,*)  ' 4 - the maximum number of main iterations is executed'
                WRITE(45,*)  ' 5 - the execution time exceeded the time limit'
                WRITE(45,*)  ' 7 - after "nochange_mainit" main iterations change in the objective is less than "nochange_tol"' 
            
        END SELECT
            
        WRITE(45,*)  ' '                
        WRITE(45,*)  'The CPU time limit: ', time_limit 
        WRITE(45,*)  ' ' 
             
        WRITE(45,*) '--------------------------------------------------------------', &
                  & '--------------------------------------------------------------', &
                  & '--------------------------------------------------------------', &
                  & '--------------------------------------------------------------'   
        WRITE(55,*) '--------------------------------------------------------------', &
                  & '--------------------------------------------------------------', &
                  & '--------------------------------------------------------------', &
                  & '--------------------------------------------------------------'                  
        
        IF (solver_ind == 2 .OR. (solver_ind == 1 .AND. optimality_condition == 2)) THEN
            WRITE(45,*) '         f1          ', 'f2          ', 'x_0      ', &
                      & '    p      ', '     m     ',   '     n  ', '|  ',  &
                      & '        n_f_tot    ', '    n_sub1_tot    ', 'n_sub2_tot ', '| ', &
                      & '      n_f_esc  ', 'n_sub_esc  ', 'esc_visits ', '| ', &
                      & '      time    ', '                      f*    ',  '                     f_0            ',  &
                      & '      termination'  
            WRITE(45,*) '                                     (number of blocks)', &
                      & '(size of block)  ' 
            WRITE(55,*) '         f1          ', 'f2          ', &
                      & '    p      ', '     m     ',   '     n  ', '|  ',  &
                      & '        n_f_ave    ', '    n_sub1_ave    ', 'n_sub2_ave ', '| ', &
                      & '      n_f_esc_ave  ', 'n_sub_esc_ave  ', 'esc_visits_ave ', '| ', &
                      & '      time_ave    ', '                      f*_ave    '  
            WRITE(55,*) '                                     (number of blocks)', &
                      & '(size of block)  ' 
        ELSE
            WRITE(45,*) '         f1          ', 'f2          ', 'x_0      ', &
                      & '    p      ', '     m     ',   '     n  ', '|  ',  &
                      & '        n_f       ', '    n_sub1       ', 'n_sub2    ', '| ', 'n_app_crit    ', '| ', &
                      & '      time    ', '                      f*    ',  '                     f_0            ',  &
                      & '      termination'     
            WRITE(45,*) '                                     (number of blocks)', &
                      & '(size of block)  ' 
            WRITE(55,*) '         f1          ', 'f2          ', &
                      & '    p      ', '     m     ',   '     n  ', '|  ',  &
                      & '        n_f_ave       ', '    n_sub1_ave       ', 'n_sub2_ave    ', '| ', &
                      & 'n_app_crit_ave    ', '| ' ,'      time_ave    ', '                      f*_ave    '
            WRITE(55,*) '                                     (number of blocks)', &
                      & '(size of block)  ' 
        END IF
             
        WRITE(45,*) '--------------------------------------------------------------', &
                  & '--------------------------------------------------------------', &
                  & '--------------------------------------------------------------', &
                  & '--------------------------------------------------------------'
        WRITE(55,*) '--------------------------------------------------------------', &
                  & '--------------------------------------------------------------', &
                  & '--------------------------------------------------------------', &
                  & '--------------------------------------------------------------'
                      
        !-------------------------------------------------------------------         
        !                           PRINT END
        !------------------------------------------------------------------- 
        
        SELECT CASE(group)
            CASE(0)
                kmax = 10
            CASE(1)
                kmax = 48
        END SELECT
        
        !-------------------------------------------------------------------         
        !                    TEST PROBLEM GROUP 0  
        !   Contains test problems 1-10 from the PBSDC paper
        !-------------------------------------------------------------------        
        IF (group == 0) THEN ! Group 0 (separable problems)

            DO i = 1, kmax

                SELECT CASE(dim_loop)   
 
                    CASE(1) ! (n= 10, 25, 50, 100, 250, 500, 1000)
                        ind_start = 1
                        ind_finish = 7

                END SELECT     

              ! Dimension of the problem          
                DO l = ind_start, ind_finish
                
                    SELECT CASE(l)
                       
                        CASE(1)
                           user_n = 10
                        
                        CASE (2)
                           user_n = 25
                        
                        CASE (3)
                           user_n = 50
                        
                        CASE (4)
                           user_n = 100
                        
                        CASE (5)
                           user_n = 250
                        
                        CASE (6)
                           user_n = 500
                        
                        CASE (7)
                           user_n = 1000

                    END SELECT  

                    ALLOCATE(problem1(user_n))
                    ALLOCATE(problem2(user_n))

                    SELECT CASE(i)   
                    
                        CASE(1)
                            problem1 = 1
                            problem2 = 2
                            bound_label1 = 0
                            bound_label2 = 0
                            
                        CASE(2)
                            problem1 = 1
                            problem2 = 3 
                            bound_label1 = 0
                            bound_label2 = 0
                            
                        CASE(3)
                            problem1 = 1
                            problem2 = 4 
                            bound_label1 = 0
                            bound_label2 = 0
                            
                        CASE(4)
                            problem1 = 2
                            problem2 = 4 
                            bound_label1 = 0
                            bound_label2 = 0
                            
                        CASE(5)
                            problem1 = 5
                            problem2 = 3 
                            bound_label1 = 0
                            bound_label2 = 0
                            
                        CASE(6)
                            problem1 = 5
                            problem2 = 4 
                            bound_label1 = 0
                            bound_label2 = 0
                            
                        CASE(7)
                            problem1 = 6
                            problem2 = 4 
                            bound_label1 = 0
                            bound_label2 = 0
                            
                        CASE(8)
                            problem1 = 7
                            problem2 = 1 
                            bound_label1 = 0
                            bound_label2 = 0
                            
                        CASE(9)
                            problem1 = 7
                            problem2 = 4 
                            bound_label1 = 0
                            bound_label2 = 0
                            
                        CASE(10)
                            problem1 = 7
                            problem2 = 8
                            bound_label1 = 0
                            bound_label2 = 0
                            
                    END SELECT   
                    
                    SELECT CASE(bound_label1)
                    
                        CASE(0)
                            ! No bounds
                        
                        CASE(1)
                            low_bound1 = -5.0_dp
                            up_bound1 = 5.0_dp
                    
                    END SELECT

                    SELECT CASE(bound_label2)
                    
                        CASE(0)
                            ! No bounds
                        
                        CASE(1)
                            low_bound2 = -5.0_dp
                            up_bound2 = 5.0_dp

                    END SELECT                  
                    
                    number_of_blocks0 = user_n  ! This holds since functions are separable (i.e. every block contains 1 variable)
                    ALLOCATE(size_of_block0(number_of_blocks0))
                    size_of_block0 = 1    ! Every block contains 1 variable
                    
                    CALL allocate_prob_data(number_of_blocks0, size_of_block0)
                    
                    f_solution_ave = 0.0_dp
                    time_ave = 0.0_dp
                    n_problems = 0
                    counter_ave = 0
                    
                    DO j = 1, 10  ! All starting points are looked through

                        IF (iprint /= 0) THEN 
                            WRITE(*,*) '------------------------------------------------------------------'
                            WRITE(*,*) '** START ** START ** START ** START ** START ** START ** START **'  
                            WRITE(*,*) '------------------------------------------------------------------'         
                            WRITE(*,*) ' '
                            WRITE(*,*) ' ', 'Problem', i, 'user_n', user_n, 'start', j
                            WRITE(*,*) ' ', 'Number of blocks', number_of_blocks0
                            WRITE(*,*) ' ', 'Size of blocks', size_of_block0(1)
                            WRITE(*,*) ' ', 'Normal start', start_normal
                            WRITE(*,*) ' '
                            WRITE(*,*) ' ', 'The first DC components of the objectives:', problem1(1)
                            WRITE(*,*) ' ', 'The Second DC components of the objectives:', problem2(1)
                            WRITE(*,*) ' '
                        END IF
                    
                        ! Correct starting points are chosen  
                        IF ((problem1(1)/=5).AND.(problem1(1)/=6).AND.(problem2(1)/=5).AND.(problem2(1)/=6)) THEN  
                   
                            start_ind = j 
                            
                            IF (solver_ind == 1) THEN  ! Method used is PBSDC
                       
                                CALL PBSDC_algorithm( f_solution, f0, mit, mrounds, &
                                            & mrounds_clarke, termination, counter, time, agg_used,  &
                                            & stepsize_used, iprint, problem1, problem2, user_n, number_of_blocks0, & 
                                            & size_of_block0, starting_points, start_ind, time_limit, &
                                            & optimality_condition, boost, nochange_mainit, nochange_tol )
                                
                            ELSE IF (solver_ind == 2) THEN   ! Method used is DBDC
                                
                                CALL DBDC_algorithm( f_solution, f0, mit, mrounds, &
                                                & mrounds_clarke, termination, counter, time, agg_used,  &
                                                & stepsize_used, iprint, problem1, problem2, user_n, & 
                                                & starting_points, start_ind, time_limit, nochange_mainit, nochange_tol )
                                
                            ELSE   ! Method used is PBDC
                            
                                CALL PBDC_algorithm( f_solution, mit, time, time_limit, &
                                                     & problem1, problem2, user_n, mrounds, & 
                                                     & termination, counter, agg_used, stepsize_used, &
                                                     & iprint, starting_points, start_ind, nochange_mainit, nochange_tol )
                            
                            END IF
                                    
                        ! Problems with exponential function has different starting points
                        ELSE  

                            start_ind = j 
                            
                            IF (solver_ind == 1) THEN  ! Method used is PBSDC
                       
                                CALL PBSDC_algorithm( f_solution, f0, mit, mrounds, &
                                        & mrounds_clarke, termination, counter, time, agg_used,  &
                                        & stepsize_used, iprint, problem1, problem2, user_n, number_of_blocks0, & 
                                        & size_of_block0, starting_points_exp, start_ind, time_limit, &
                                        & optimality_condition, boost, nochange_mainit, nochange_tol )
                                
                            ELSE IF (solver_ind == 2) THEN   ! Method used is DBDC
                                
                                CALL DBDC_algorithm( f_solution, f0, mit, mrounds, &
                                        & mrounds_clarke, termination, counter, time, agg_used,  &
                                        & stepsize_used, iprint, problem1, problem2, user_n, & 
                                        & starting_points_exp, start_ind, time_limit, nochange_mainit, nochange_tol )
                                
                            ELSE   ! Method used is PBDC
                            
                                CALL PBDC_algorithm( f_solution, mit, time, time_limit, &
                                                     & problem1, problem2, user_n, mrounds, & 
                                                     & termination, counter, agg_used, stepsize_used, &
                                                     & iprint, starting_points_exp, start_ind, nochange_mainit, nochange_tol )
                            
                            END IF

                        END IF  

                        IF (iprint /= 0) THEN 
                            WRITE(*,*) ' '
                            WRITE(*,*) '------------------------------------------------------------------'
                            WRITE(*,*) '** END ** END ** END ** END ** END ** END ** END ** END ** END **'  
                            WRITE(*,*) '------------------------------------------------------------------'
                            WRITE(*,*) ' '      
                        END IF
                   
                        IF (solver_ind == 2 .OR. (solver_ind == 1 .AND. optimality_condition == 2)) THEN
                            WRITE(45,*) problem1(1), problem2(1), start_ind, number_of_blocks0, size_of_block0(1), user_n, &
                                        & '| ', counter(3)+counter(7), counter(4)+counter(8), counter(5)+counter(8),  '| ', &
                                        & counter(7), counter(8), counter(6), '| ',  & 
                                        & time, f_solution, f0, termination
                        ELSE 
                            WRITE(45,*) problem1(1), problem2(1), start_ind, number_of_blocks0, size_of_block0(1), user_n, &
                                        & '| ', counter(3), counter(4), counter(5),  '| ', counter(6), '| ', &
                                        & time, f_solution, f0, termination
                        END IF
                        
                        f_solution_ave = f_solution_ave + f_solution
                        time_ave = time_ave + time
                        n_problems = n_problems + 1
                        counter_ave = counter_ave + counter
                     
                    END DO   

                    c = 1.0_dp / n_problems                 

                    IF (solver_ind == 2 .OR. (solver_ind == 1 .AND. optimality_condition == 2)) THEN
                        WRITE(55,*) problem1(1), problem2(1), number_of_blocks0, size_of_block0(1), user_n, &
                                    & '| ', c*(counter_ave(3)+counter_ave(7)), c*(counter_ave(4)+counter_ave(8)), &
                                    & c*(counter_ave(5)+counter_ave(8)),  '| ', &
                                    & c*counter_ave(7), c*counter_ave(8), c*counter_ave(6), '| ',  & 
                                    & c*time_ave, c*f_solution_ave
                    ELSE
                        WRITE(55,*) problem1(1), problem2(1), number_of_blocks0, size_of_block0(1), user_n, &
                                    & '| ', c*counter_ave(3), c*counter_ave(4), c*counter_ave(5),  '| ', &
                                    & c*counter_ave(6), '| ', c*time_ave, c*f_solution_ave
                    END IF                  
                    
                    !-----------------------------------
                    ! Deallocations
                    !-----------------------------------
                    
                    CALL deallocate_prob_data()
                    DEALLOCATE(size_of_block0)
                    DEALLOCATE(problem1)
                    DEALLOCATE(problem2)
                    
                END DO
        
            END DO  

        ELSE ! Group = 1 (block-separable problems)
        
            DO i = 1, kmax

                SELECT CASE(block_number_label)  
         
                    CASE(1) ! (block numbers = 10, 30, 50)
                        ind_start_block_number = 1
                        ind_finish_block_number = 3

                    CASE(2) ! (block numbers = 80, 100, 120)
                        ind_start_block_number = 4
                        ind_finish_block_number = 6

                END SELECT
                    
                DO l = ind_start_block_number, ind_finish_block_number
                        
                    SELECT CASE(l)
                           
                        CASE(1)
                            number_of_blocks0 = 10
                            
                        CASE(2)
                            number_of_blocks0 = 30
                            
                        CASE(3)
                            number_of_blocks0 = 50
                            
                        CASE(4)
                            number_of_blocks0 = 80
                            
                        CASE(5)
                            number_of_blocks0 = 100
                            
                        CASE(6)
                            number_of_blocks0 = 120
                            
                    END SELECT
                    
                    ALLOCATE(size_of_block0(number_of_blocks0))
                    ALLOCATE(problem1(number_of_blocks0))
                    ALLOCATE(problem2(number_of_blocks0))
                    
                    SELECT CASE(group)  
                
                        CASE(1) ! ---- Group 1 of Test Problems solved -----
                        !-------------------------------------------------------------------         
                        !                    TEST PROBLEM GROUP 1  
                        !   Contains test problems 11-58 from the PBSDC paper
                        !------------------------------------------------------------------- 
                         
                            SELECT CASE(i)   
             
                                CASE(1)
                                    problem1 = 1001
                                    problem2 = 1003  
                                    
                                CASE(2)
                                    problem1 = 1001
                                    problem2 = 1005  
                                    
                                CASE(3)
                                    problem1 = 1001
                                    problem2 = 1006  

                                CASE(4)
                                    problem1 = 1001
                                    problem2 = 1007  

                                CASE(5)
                                    problem1 = 1001
                                    problem2 = 1012

                                CASE(6)
                                    problem1 = 1002
                                    problem2 = 1005

                                CASE(7)
                                    problem1 = 1002
                                    problem2 = 1007

                                CASE(8)
                                    problem1 = 1003
                                    problem2 = 1002

                                CASE(9)
                                    problem1 = 1003
                                    problem2 = 1004

                                CASE(10)
                                    problem1 = 1003
                                    problem2 = 1005

                                CASE(11)
                                    problem1 = 1003
                                    problem2 = 1006

                                CASE(12)
                                    problem1 = 1004
                                    problem2 = 1005

                                CASE(13)
                                    problem1 = 1004
                                    problem2 = 1007

                                CASE(14)
                                    problem1 = 1005
                                    problem2 = 1007

                                CASE(15)
                                    problem1 = 1008
                                    problem2 = 1006

                                CASE(16)
                                    problem1 = 1008
                                    problem2 = 1007

                                CASE(17)
                                    problem1 = 1008
                                    problem2 = 1013

                                CASE(18)
                                    problem1 = 1009
                                    problem2 = 1007

                                CASE(19)
                                    problem1 = 1009
                                    problem2 = 1012

                                CASE(20)
                                    problem1 = 1009
                                    problem2 = 1013

                                CASE(21)
                                    problem1 = 1010
                                    problem2 = 1007

                                CASE(22)
                                    problem1 = 1010
                                    problem2 = 1012

                                CASE(23)
                                    problem1 = 1010
                                    problem2 = 1013

                                CASE(24)
                                    problem1 = 1011
                                    problem2 = 1003

                                CASE(25)
                                    problem1 = 1011
                                    problem2 = 1006

                                CASE(26)
                                    problem1 = 1011
                                    problem2 = 1008

                                CASE(27)
                                    problem1 = 1011
                                    problem2 = 1012

                                CASE(28)
                                    problem1 = 1011
                                    problem2 = 1013

                                CASE(29)
                                    problem1 = 1011
                                    problem2 = 1015

                                CASE(30)
                                    problem1 = 1011
                                    problem2 = 1016

                                CASE(31)
                                    problem1 = 1014
                                    problem2 = 1001

                                CASE(32)
                                    problem1 = 1014
                                    problem2 = 1002

                                CASE(33)
                                    problem1 = 1014
                                    problem2 = 1003

                                CASE(34)
                                    problem1 = 1014
                                    problem2 = 1004

                                CASE(35)
                                    problem1 = 1014
                                    problem2 = 1005

                                CASE(36)
                                    problem1 = 1014
                                    problem2 = 1006

                                CASE(37)
                                    problem1 = 1014
                                    problem2 = 1007

                                CASE(38)
                                    problem1 = 1014
                                    problem2 = 1008

                                CASE(39)
                                    problem1 = 1014
                                    problem2 = 1009

                                CASE(40)
                                    problem1 = 1014
                                    problem2 = 1012

                                CASE(41)
                                    problem1 = 1014
                                    problem2 = 1013

                                CASE(42)
                                    problem1 = 1014
                                    problem2 = 1016

                                CASE(43)
                                    problem1 = 1015
                                    problem2 = 1002

                                CASE(44)
                                    problem1 = 1015
                                    problem2 = 1004

                                CASE(45)
                                    problem1 = 1015
                                    problem2 = 1006

                                CASE(46)
                                    problem1 = 1015
                                    problem2 = 1013

                                CASE(47)
                                    problem1 = 1015
                                    problem2 = 1016

                                CASE(48)
                                    problem1 = 1016
                                    problem2 = 1006                                     
                                    
                            END SELECT  
                            
                    END SELECT
                    
                    SELECT CASE(block_size_label)   
     
                        CASE(1) ! (block size = 5, 10, 20)
                            ind_start_block_size = 1
                            ind_finish_block_size = 3

                        CASE(2) ! (block size = 25, 50, 100)
                            ind_start_block_size = 4
                            ind_finish_block_size = 6

                    END SELECT     

                    ! Block size of the problem          
                    DO k = ind_start_block_size, ind_finish_block_size
                    
                        SELECT CASE(k)
                           
                            CASE(1)
                               size_of_block0 = 5
                            
                            CASE(2)
                               size_of_block0 = 10
                            
                            CASE(3)
                               size_of_block0 = 20
                            
                            CASE(4)
                               size_of_block0 = 25
                            
                            CASE(5)
                               size_of_block0 = 50
                            
                            CASE(6)
                               size_of_block0 = 100
                            
                        END SELECT  

                        user_n = 0
                        DO j = 1,number_of_blocks0
                            user_n = user_n + size_of_block0(j)
                        END DO
                        
                        CALL allocate_prob_data(number_of_blocks0, size_of_block0)
                        
                        f_solution_ave = 0.0_dp
                        time_ave = 0.0_dp
                        n_problems = 0
                        counter_ave = 0

                        DO j = 1, 10  ! All starting points are looked through
                        
                            SELECT CASE(group)
                                CASE(1)
                                    start_normal = ( (problem1(1) .NE. 1009) .AND. (problem1(1) .NE. 1010) &
                                                     & .AND. (problem1(1) .NE. 1014) .AND. (problem2(1) .NE. 1009) &
                                                     & .AND. (problem2(1) .NE. 1010) .AND. (problem2(1) .NE. 1014) )
                            END SELECT
                            
                            IF (iprint /= 0) THEN 
                                WRITE(*,*) '------------------------------------------------------------------'
                                WRITE(*,*) '** START ** START ** START ** START ** START ** START ** START **'  
                                WRITE(*,*) '------------------------------------------------------------------'         
                                WRITE(*,*) ' '
                                WRITE(*,*) ' ', 'Problem', i, 'user_n', user_n, 'start', j
                                WRITE(*,*) ' ', 'Number of blocks', number_of_blocks0
                                WRITE(*,*) ' ', 'Size of blocks', size_of_block0(1)
                                WRITE(*,*) ' ', 'Normal start', start_normal
                                WRITE(*,*) ' '
                                WRITE(*,*) ' ', 'The first DC components of the objectives:', problem1(1)
                                WRITE(*,*) ' ', 'The Second DC components of the objectives:', problem2(1)
                                WRITE(*,*) ' '
                            END IF

                            ! Correct starting points are chosen  
                            IF (start_normal) THEN  
                            
                                start_ind = j
                            
                                IF (solver_ind == 1) THEN  ! Method used is PBSDC
                                    
                                    IF (iprint >= 1) THEN         
                                        WRITE(*,*) ' ', 'boosting', boost, '1=no, 2=yes'
                                        WRITE(*,*) ' '
                                    END IF
                       
                                    CALL PBSDC_algorithm( f_solution, f0, mit, mrounds, &
                                                & mrounds_clarke, termination, counter, time, agg_used,  &
                                                & stepsize_used, iprint, problem1, problem2, user_n, number_of_blocks0, & 
                                                & size_of_block0, starting_points, start_ind, time_limit, &
                                                & optimality_condition, boost, nochange_mainit, nochange_tol )
                                
                                ELSE IF (solver_ind == 2) THEN   ! Method used is DBDC
                                    
                                    CALL DBDC_algorithm( f_solution, f0, mit, mrounds, &
                                                    & mrounds_clarke, termination, counter, time, agg_used,  &
                                                    & stepsize_used, iprint, problem1, problem2, user_n, & 
                                                    & starting_points, start_ind, time_limit, nochange_mainit, nochange_tol )
                                    
                                ELSE   ! Method used is PBDC
                                
                                    CALL PBDC_algorithm( f_solution, mit, time, time_limit, &
                                                         & problem1, problem2, user_n, mrounds, & 
                                                         & termination, counter, agg_used, stepsize_used, &
                                                         & iprint, starting_points, start_ind, nochange_mainit, nochange_tol )
                            
                                END IF
                                       
                            ELSE  ! Different starting points due to exponential function
                                
                                start_ind = j
                            
                                IF (solver_ind == 1) THEN  ! Method used is PBSDC
                                
                                    IF (iprint >= 1) THEN         
                                        WRITE(*,*) ' ', 'boosting', boost, '1=no, 2=yes'
                                        WRITE(*,*) ' '
                                    END IF
                       
                                    CALL PBSDC_algorithm( f_solution, f0, mit, mrounds, &
                                                & mrounds_clarke, termination, counter, time, agg_used,  &
                                                & stepsize_used, iprint, problem1, problem2, user_n, number_of_blocks0, & 
                                                & size_of_block0, starting_points_exp, start_ind, time_limit, &
                                                & optimality_condition, boost, nochange_mainit, nochange_tol )
                                
                                ELSE IF (solver_ind == 2) THEN   ! Method used is DBDC
                                    
                                    CALL DBDC_algorithm( f_solution, f0, mit, mrounds, &
                                                    & mrounds_clarke, termination, counter, time, agg_used,  &
                                                    & stepsize_used, iprint, problem1, problem2, user_n, & 
                                                    & starting_points_exp, start_ind, time_limit, nochange_mainit, nochange_tol )
                                    
                                ELSE   ! Method used is PBDC
                                
                                    CALL PBDC_algorithm( f_solution, mit, time, time_limit, &
                                                         & problem1, problem2, user_n, mrounds, & 
                                                         & termination, counter, agg_used, stepsize_used, &
                                                         & iprint, starting_points_exp, start_ind, nochange_mainit, nochange_tol )
                            
                                END IF

                            END IF   

                            IF (iprint /= 0) THEN 
                                WRITE(*,*) ' '
                                WRITE(*,*) '------------------------------------------------------------------'
                                WRITE(*,*) '** END ** END ** END ** END ** END ** END ** END ** END ** END **'  
                                WRITE(*,*) '------------------------------------------------------------------'
                                WRITE(*,*) ' '      
                            END IF
                            
                            IF (solver_ind == 2 .OR. (solver_ind == 1 .AND. optimality_condition == 2)) THEN
                                WRITE(45,*) problem1(1), problem2(1), start_ind, number_of_blocks0, size_of_block0(1), user_n, &
                                         & '| ', counter(3)+counter(7), counter(4)+counter(8), counter(5)+counter(8),  '| ', &
                                         & counter(7), counter(8), counter(6), '| ',  & 
                                         & time, f_solution, f0, termination
                            ELSE 
                                WRITE(45,*) problem1(1), problem2(1), start_ind, number_of_blocks0, size_of_block0(1), user_n, &
                                         & '| ', counter(3), counter(4), counter(5),  '| ', counter(6), '| ', &
                                         & time, f_solution, f0, termination
                            END IF

                            f_solution_ave = f_solution_ave + f_solution
                            time_ave = time_ave + time
                            n_problems = n_problems + 1
                            counter_ave = counter_ave + counter
                         
                        END DO   

                        c = 1.0_dp / n_problems                 

                        IF (solver_ind == 2 .OR. (solver_ind == 1 .AND. optimality_condition == 2)) THEN
                            WRITE(55,*) problem1(1), problem2(1), number_of_blocks0, size_of_block0(1), user_n, &
                                        & '| ', c*(counter_ave(3)+counter_ave(7)), c*(counter_ave(4)+counter_ave(8)), &
                                        & c*(counter_ave(5)+counter_ave(8)),  '| ', &
                                        & c*counter_ave(7), c*counter_ave(8), c*counter_ave(6), '| ',  & 
                                        & c*time_ave, c*f_solution_ave
                        ELSE
                            WRITE(55,*) problem1(1), problem2(1), number_of_blocks0, size_of_block0(1), user_n, &
                                        & '| ', c*counter_ave(3), c*counter_ave(4), c*counter_ave(5),  '| ', &
                                        & c*counter_ave(6), '| ', c*time_ave, c*f_solution_ave
                        END IF
                        
                        CALL deallocate_prob_data()
                        
                    END DO  

                    DEALLOCATE(size_of_block0)
                    DEALLOCATE(problem1)
                    DEALLOCATE(problem2)
                    
                END DO    

            END DO  
            
        END IF          
        
        CLOSE(45)
        CLOSE(55)
             
      END PROGRAM tpbsdc