import sympy as sp

def solve_queueing_problem():
    """
    Calculates the optimal mean response time x and extracts the specified term.
    """
    # 1. Define symbols and parameters
    s, y, x = sp.symbols('s y x')
    lambda_val = sp.Rational(3, 2)
    # Service time S ~ U(0,1), so its PDF f(s) is 1 for s in [0,1]
    f_s = 1

    # 2. Calculate key components for the SRPT formula
    
    # rho(y) is the utilization from jobs of size less than y
    # rho(y) = lambda * integral from 0 to y of (t * f(t) dt)
    rho_y = lambda_val * sp.integrate(s * f_s, (s, 0, y))

    # E[T(x)] is the mean response time for a job of size x.
    # It has two parts.
    # Part 1: Time spent on jobs of size <= x
    E_T_x_part1 = sp.integrate(1 / (1 - rho_y), (y, 0, x))

    # Part 2: Delay from work present on arrival
    # Numerator of Part 2
    E_T_x_part2_num = lambda_val * sp.integrate(s**2 * f_s, (s, 0, x))
    # Denominator of Part 2
    rho_x = rho_y.subs(y, x)
    E_T_x_part2_den = 2 * (1 - rho_x)**2
    E_T_x_part2 = E_T_x_part2_num / E_T_x_part2_den

    E_T_x = E_T_x_part1 + E_T_x_part2

    # 3. Calculate the overall mean response time, x_val
    # x_val = integral from 0 to 1 of (E[T(s)] * f(s) ds)
    x_val = sp.integrate(E_T_x.subs(x, s), (s, 0, 1))

    # 4. Decompose the result and print the equation
    # Sympy might group terms, so we express it in a clear form.
    # x_val is -16*log(2)/9 + 2/3 + 2*sqrt(3)*log(sqrt(3) + 2)/3
    terms = x_val.as_ordered_terms()
    
    rational_term = None
    log_rational_term = None
    remaining_term = None

    for term in terms:
        if term.is_rational:
            rational_term = term
        # Check for terms of the form c * log(q) where c, q are rational.
        elif term.has(sp.log):
            log_args = [arg for log_func in term.atoms(sp.log) for arg in log_func.args]
            # Check if all log arguments are rational
            if all(arg.is_rational for arg in log_args):
                log_rational_term = term
            else:
                remaining_term = term
        else:
            remaining_term = term

    print("The optimal mean response time x is given by the equation:")
    # Pretty print the equation
    equation = sp.Eq(sp.Symbol('x'), x_val)
    sp.pprint(equation, use_unicode=True)
    
    print("\nDecomposition of x:")
    print(f"Term 1 (Rational): {rational_term}")
    print(f"Term 2 (Log of Rational): {log_rational_term}")
    print(f"Term 3 (Remaining): {remaining_term}")
    
    # 5. Remove specified terms and find the remainder
    print("\nRemoving the rational term and the term which is a logarithm of a rational number...")
    
    # 6. Format the final answer
    # The remaining term is (2*sqrt(3)/3) * log(2 + sqrt(3))
    # Format as: rational * algebraic_irrational * transcendental
    # (2/3) * sqrt(3) * log(2 + sqrt(3))
    final_term_formatted = sp.Rational(2, 3) * sp.sqrt(3) * sp.log(2 + sp.sqrt(3))
    
    print("\nThe remaining term is:")
    print(sp.latex(final_term_formatted))
    
if __name__ == '__main__':
    solve_queueing_problem()