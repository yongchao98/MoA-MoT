import sympy
from sympy import Symbol, integrate, log, sqrt, Rational, latex, Add, Mul

def solve_queueing_problem():
    """
    Solves the queueing theory problem symbolically to find the specific term
    of the optimal mean response time x.
    This script calculates the mean response time for an M/G/1 queue with SRPT scheduling,
    then isolates the required term based on the problem description.
    """
    # Define symbols and parameters for the M/G/1 SRPT queue.
    s = Symbol('s')
    y = Symbol('y')
    z = Symbol('z')  # Dummy variable for integration

    lambda_val = Rational(3, 2)
    # For U(0, 1), the PDF f_S(s) is 1 on the interval [0, 1].
    f_s = 1

    # The mean response time `x` is the integral of E[T(s)] over [0, 1].
    # E[T(s)] = E[W(s)] + s
    # We calculate E[W(s)] first, which is the sum of two terms.

    # rho(s): traffic intensity from jobs of size <= s
    rho_s = lambda_val * integrate(y * f_s, (y, 0, s))

    # First term of E[W(s)]
    term1_num = lambda_val * integrate(y**2 * f_s, (y, 0, s))
    term1_den = 2 * (1 - rho_s)**2
    term1_s = term1_num / term1_den

    # Second term of E[W(s)]
    # We need to define rho(y) inside this integral
    rho_y = lambda_val * integrate(z * f_s, (z, 0, y))
    integrand_term2 = (lambda_val * y * f_s) / (1 - rho_y)
    term2_s = integrate(integrand_term2, (y, 0, s))

    # E[T(s)]: mean response time for a job of size s
    E_T_s = term1_s + term2_s + s

    # Calculate x = E[T], the overall mean response time
    # The integration is over the support of the job size distribution, [0, 1].
    x = integrate(E_T_s * f_s, (s, 0, 1))

    # Sympy might return expressions with atanh; we convert to log for clarity.
    x_expanded = x.expand().rewrite(log)
    
    # Identify the components of x as per the problem statement
    if isinstance(x_expanded, Add):
        all_terms = x_expanded.args
    else:
        all_terms = [x_expanded]
        
    rational_terms = []
    log_rational_terms = []
    remaining_terms = []
    
    for term in all_terms:
        if term.is_rational:
            rational_terms.append(term)
            continue
        
        # Heuristic to check for terms of the form q * log(r)
        # where q and r are rational.
        is_log_rational = False
        if term.is_Mul:
            coeff = Mul(*[arg for arg in term.args if not isinstance(arg, log)])
            log_part = [arg for arg in term.args if isinstance(arg, log)]
            if len(log_part) == 1 and coeff.is_rational:
                log_arg = log_part[0].args[0]
                if log_arg.is_rational:
                    is_log_rational = True

        if is_log_rational:
            log_rational_terms.append(term)
        else:
            remaining_terms.append(term)

    # Combine the parts for printing
    rational_part = Add(*rational_terms)
    log_rational_part = Add(*log_rational_terms)
    remaining_part = Add(*remaining_terms)
    
    # Output the required numbers and the final answer
    print("The optimal mean response time `x` is composed of the following parts:")
    # Use replace to show ln instead of log for natural logarithm
    print(f"  x = ({latex(rational_part)}) + ({latex(log_rational_part).replace(r'log', r'ln')}) + ({latex(remaining_part).replace(r'log', r'ln')})")
    
    print("\nBased on the problem, we remove the following terms:")
    print(f"1. The additive rational term: {latex(rational_part)}")
    print(f"2. The additive term which is a logarithm of a rational number: {latex(log_rational_part).replace(r'log', r'ln')}")

    print("\nThe remaining term is:")
    # Format the final answer in LaTeX
    final_latex_str = latex(remaining_part).replace(r'log', r'ln')
    print(final_latex_str)

if __name__ == '__main__':
    solve_queueing_problem()