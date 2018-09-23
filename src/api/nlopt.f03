module nlopt
 use, intrinsic :: iso_c_binding
 implicit none
 public

 enum, bind(c)
  ! FIXME this list should be generated from nlopt.h
  enumerator :: NLOPT_GN_DIRECT = 0, &
    NLOPT_GN_DIRECT_L, &
    NLOPT_GN_DIRECT_L_RAND, &
    NLOPT_GN_DIRECT_NOSCAL, &
    NLOPT_GN_DIRECT_L_NOSCAL, &
    NLOPT_GN_DIRECT_L_RAND_NOSCAL, &
    NLOPT_GN_ORIG_DIRECT, &
    NLOPT_GN_ORIG_DIRECT_L, &
    NLOPT_GD_STOGO, &
    NLOPT_GD_STOGO_RAND, &
    NLOPT_LD_LBFGS_NOCEDAL, &
    NLOPT_LD_LBFGS, &
    NLOPT_LN_PRAXIS, &
    NLOPT_LD_VAR1, &
    NLOPT_LD_VAR2, &
    NLOPT_LD_TNEWTON, &
    NLOPT_LD_TNEWTON_RESTART, &
    NLOPT_LD_TNEWTON_PRECOND, &
    NLOPT_LD_TNEWTON_PRECOND_RESTART, &
    NLOPT_GN_CRS2_LM, &
    NLOPT_GN_MLSL, &
    NLOPT_GD_MLSL, &
    NLOPT_GN_MLSL_LDS, &
    NLOPT_GD_MLSL_LDS, &
    NLOPT_LD_MMA, &
    NLOPT_LN_COBYLA, &
    NLOPT_LN_NEWUOA, &
    NLOPT_LN_NEWUOA_BOUND, &
    NLOPT_LN_NELDERMEAD, &
    NLOPT_LN_SBPLX, &
    NLOPT_LN_AUGLAG, &
    NLOPT_LD_AUGLAG, &
    NLOPT_LN_AUGLAG_EQ, &
    NLOPT_LD_AUGLAG_EQ, &
    NLOPT_LN_BOBYQA, &
    NLOPT_GN_ISRES, &
    NLOPT_AUGLAG, &
    NLOPT_AUGLAG_EQ, &
    NLOPT_G_MLSL, &
    NLOPT_G_MLSL_LDS, &
    NLOPT_LD_SLSQP, &
    NLOPT_LD_CCSAQ, &
    NLOPT_GN_ESCH, &
    NLOPT_GN_AGS, &
    NLOPT_NUM_ALGORITHMS
   enumerator ::  NLOPT_FAILURE = -1, &
    NLOPT_INVALID_ARGS = -2, &
    NLOPT_OUT_OF_MEMORY = -3, &
    NLOPT_ROUNDOFF_LIMITED = -4, &
    NLOPT_FORCED_STOP = -5, &
    NLOPT_SUCCESS = 1, &
    NLOPT_MINF_MAX_REACHED = 2, &
    NLOPT_STOPVAL_REACHED = 2, &
    NLOPT_FTOL_REACHED = 3, &
    NLOPT_XTOL_REACHED = 4, &
    NLOPT_MAXEVAL_REACHED = 5, &
    NLOPT_MAXTIME_REACHED = 6
 end enum

 abstract interface
 !typedef double (*nlopt_func) (unsigned n, const double *x,
 !                              double *gradient, /* NULL if not needed */
 !                              void *func_data);
 function nlopt_func(n, x, gradient, func_data) result(ret) bind(c)
  import c_double, c_ptr, c_int
  real(c_double) :: ret, x(*), gradient(*)
  integer(c_int), value :: n
  type(c_ptr) :: func_data
 end function

 !typedef void (*nlopt_mfunc) (unsigned m, double *result, unsigned n, const double *x,
 !                             double *gradient, /* NULL if not needed */
 !                             void *func_data);
 subroutine nlopt_mfunc(m, result, n, x, gradient, func_data) bind(c)
  import c_double, c_int, c_ptr
  real(c_double) :: result(*), x(*), gradient(*)
  integer(c_int), value :: m, n
  type(c_ptr) :: func_data
 end subroutine

 !typedef void (*nlopt_precond) (unsigned n, const double *x, const double *v, double *vpre, void *data);
 subroutine nlopt_precond(n, x, v, vpre, data) bind(c)
  import c_double, c_int, c_ptr
  real(c_double) :: x(*), v(*), vpre(*)
  integer(c_int), value :: n
  type(c_ptr) :: data
 end subroutine
 end interface

 !struct nlopt_opt_s;             /* opaque structure, defined internally */
 !typedef struct nlopt_opt_s *nlopt_opt;
 type, bind(c) :: nlopt_opt
  private
  type(c_ptr) :: ptr
 end type
 ! FIXME: is this the best way to describe a typed opaque pointer? is it guaranteed to be compatible?

 interface

 !NLOPT_EXTERN(const char *) nlopt_algorithm_name(nlopt_algorithm a);
 function nlopt_algorithm_name(a) result(name) bind(c)
  import c_ptr
  integer(kind(NLOPT_NUM_ALGORITHMS)), value :: a ! FIXME: may not be C interoperable?
  type(c_ptr) :: name
  ! maybe it makes sense to write a wrapper that would make a character,pointer,dimension(:) :: string by calling call c_f_pointer(ptr, string)
 end function

 !NLOPT_EXTERN(void) nlopt_srand(unsigned long seed);
 subroutine nlopt_stand(seed) bind(c)
  import c_long
  integer(c_long), value :: seed
 end subroutine

 !NLOPT_EXTERN(void) nlopt_srand_time(void);
 subroutine nlopt_srand_time() bind(c)
 end subroutine

 !NLOPT_EXTERN(void) nlopt_version(int *major, int *minor, int *bugfix);
 subroutine nlopt_version(major, minor, bugfix) bind(c)
  import c_int
  integer(c_int), intent(out) :: major, minor, bugfix
 end subroutine

!/*************************** OBJECT-ORIENTED API **************************/
!/* The style here is that we create an nlopt_opt "object" (an opaque pointer),
!   then set various optimization parameters, and then execute the
!   algorithm.  In this way, we can add more and more optimization parameters
!   (including algorithm-specific ones) without breaking backwards
!   compatibility, having functions with zillions of parameters, or
!   relying non-reentrantly on global variables.*/

!/* the only immutable parameters of an optimization are the algorithm and
!   the dimension n of the problem, since changing either of these could
!   have side-effects on lots of other parameters */
 !NLOPT_EXTERN(nlopt_opt) nlopt_create(nlopt_algorithm algorithm, unsigned n);
 function nlopt_create(algo, n) bind(c) result(ret)
  import nlopt_opt, c_int
  type(nlopt_opt) :: ret
  integer(kind(NLOPT_NUM_ALGORITHMS)), value :: algo
  integer(c_int), value :: n
 end function

 !NLOPT_EXTERN(void) nlopt_destroy(nlopt_opt opt);
 subroutine nlopt_destroy(opt) bind(c)
  import nlopt_opt
  type(nlopt_opt), value :: opt
 end subroutine
!NLOPT_EXTERN(nlopt_opt) nlopt_copy(const nlopt_opt opt);
!
!NLOPT_EXTERN(nlopt_result) nlopt_optimize(nlopt_opt opt, double *x, double *opt_f);
!
!NLOPT_EXTERN(nlopt_result) nlopt_set_min_objective(nlopt_opt opt, nlopt_func f, void *f_data);
!NLOPT_EXTERN(nlopt_result) nlopt_set_max_objective(nlopt_opt opt, nlopt_func f, void *f_data);
!
!NLOPT_EXTERN(nlopt_result) nlopt_set_precond_min_objective(nlopt_opt opt, nlopt_func f, nlopt_precond pre, void *f_data);
!NLOPT_EXTERN(nlopt_result) nlopt_set_precond_max_objective(nlopt_opt opt, nlopt_func f, nlopt_precond pre, void *f_data);
!
!NLOPT_EXTERN(nlopt_algorithm) nlopt_get_algorithm(const nlopt_opt opt);
!NLOPT_EXTERN(unsigned) nlopt_get_dimension(const nlopt_opt opt);
!
!NLOPT_EXTERN(const char *) nlopt_get_errmsg(nlopt_opt opt);
!
!/* constraints: */
!
!NLOPT_EXTERN(nlopt_result) nlopt_set_lower_bounds(nlopt_opt opt, const double *lb);
!NLOPT_EXTERN(nlopt_result) nlopt_set_lower_bounds1(nlopt_opt opt, double lb);
!NLOPT_EXTERN(nlopt_result) nlopt_get_lower_bounds(const nlopt_opt opt, double *lb);
!NLOPT_EXTERN(nlopt_result) nlopt_set_upper_bounds(nlopt_opt opt, const double *ub);
!NLOPT_EXTERN(nlopt_result) nlopt_set_upper_bounds1(nlopt_opt opt, double ub);
!NLOPT_EXTERN(nlopt_result) nlopt_get_upper_bounds(const nlopt_opt opt, double *ub);
!
!NLOPT_EXTERN(nlopt_result) nlopt_remove_inequality_constraints(nlopt_opt opt);
!NLOPT_EXTERN(nlopt_result) nlopt_add_inequality_constraint(nlopt_opt opt, nlopt_func fc, void *fc_data, double tol);
!NLOPT_EXTERN(nlopt_result) nlopt_add_precond_inequality_constraint(nlopt_opt opt, nlopt_func fc, nlopt_precond pre, void *fc_data, double tol);
!NLOPT_EXTERN(nlopt_result) nlopt_add_inequality_mconstraint(nlopt_opt opt, unsigned m, nlopt_mfunc fc, void *fc_data, const double *tol);
!
!NLOPT_EXTERN(nlopt_result) nlopt_remove_equality_constraints(nlopt_opt opt);
!NLOPT_EXTERN(nlopt_result) nlopt_add_equality_constraint(nlopt_opt opt, nlopt_func h, void *h_data, double tol);
!NLOPT_EXTERN(nlopt_result) nlopt_add_precond_equality_constraint(nlopt_opt opt, nlopt_func h, nlopt_precond pre, void *h_data, double tol);
!NLOPT_EXTERN(nlopt_result) nlopt_add_equality_mconstraint(nlopt_opt opt, unsigned m, nlopt_mfunc h, void *h_data, const double *tol);
!
!/* stopping criteria: */
!
!NLOPT_EXTERN(nlopt_result) nlopt_set_stopval(nlopt_opt opt, double stopval);
!NLOPT_EXTERN(double) nlopt_get_stopval(const nlopt_opt opt);
!
!NLOPT_EXTERN(nlopt_result) nlopt_set_ftol_rel(nlopt_opt opt, double tol);
!NLOPT_EXTERN(double) nlopt_get_ftol_rel(const nlopt_opt opt);
!NLOPT_EXTERN(nlopt_result) nlopt_set_ftol_abs(nlopt_opt opt, double tol);
!NLOPT_EXTERN(double) nlopt_get_ftol_abs(const nlopt_opt opt);
!
!NLOPT_EXTERN(nlopt_result) nlopt_set_xtol_rel(nlopt_opt opt, double tol);
!NLOPT_EXTERN(double) nlopt_get_xtol_rel(const nlopt_opt opt);
!NLOPT_EXTERN(nlopt_result) nlopt_set_xtol_abs1(nlopt_opt opt, double tol);
!NLOPT_EXTERN(nlopt_result) nlopt_set_xtol_abs(nlopt_opt opt, const double *tol);
!NLOPT_EXTERN(nlopt_result) nlopt_get_xtol_abs(const nlopt_opt opt, double *tol);
!
!NLOPT_EXTERN(nlopt_result) nlopt_set_maxeval(nlopt_opt opt, int maxeval);
!NLOPT_EXTERN(int) nlopt_get_maxeval(const nlopt_opt opt);
!
!NLOPT_EXTERN(int) nlopt_get_numevals(const nlopt_opt opt);
!
!NLOPT_EXTERN(nlopt_result) nlopt_set_maxtime(nlopt_opt opt, double maxtime);
!NLOPT_EXTERN(double) nlopt_get_maxtime(const nlopt_opt opt);
!
!NLOPT_EXTERN(nlopt_result) nlopt_force_stop(nlopt_opt opt);
!NLOPT_EXTERN(nlopt_result) nlopt_set_force_stop(nlopt_opt opt, int val);
!NLOPT_EXTERN(int) nlopt_get_force_stop(const nlopt_opt opt);
!
!/* more algorithm-specific parameters */
!
!NLOPT_EXTERN(nlopt_result) nlopt_set_local_optimizer(nlopt_opt opt, const nlopt_opt local_opt);
!
!NLOPT_EXTERN(nlopt_result) nlopt_set_population(nlopt_opt opt, unsigned pop);
!NLOPT_EXTERN(unsigned) nlopt_get_population(const nlopt_opt opt);
!
!NLOPT_EXTERN(nlopt_result) nlopt_set_vector_storage(nlopt_opt opt, unsigned dim);
!NLOPT_EXTERN(unsigned) nlopt_get_vector_storage(const nlopt_opt opt);
!
!NLOPT_EXTERN(nlopt_result) nlopt_set_default_initial_step(nlopt_opt opt, const double *x);
!NLOPT_EXTERN(nlopt_result) nlopt_set_initial_step(nlopt_opt opt, const double *dx);
!NLOPT_EXTERN(nlopt_result) nlopt_set_initial_step1(nlopt_opt opt, double dx);
!NLOPT_EXTERN(nlopt_result) nlopt_get_initial_step(const nlopt_opt opt, const double *x, double *dx);
 end interface
end module
