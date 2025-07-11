import sympy

def solve_queueing_problem():
    """
    This function calculates the optimal mean response time 'x' for the given
    M/G/1 SRPT queue and extracts the specified remaining term.
    """
    # Define symbols and constants
    s, y = sympy.symbols('s y')
    lambda_rate = sympy.Rational(3, 2)

    # Job sizes S are from a Uniform distribution U(0, 1).
    # The PDF f(s) is 1 for 0 <= s <= 1.
    # The CDF F(s) is s for 0 <= s <= 1.
    # dF(s) is ds.

    # 1. Calculate Mean Service Time E[S]
    E_S = sympy.integrate(s, (s, 0, 1))

    # 2. Define components for the mean waiting time calculation
    # rho_s: traffic intensity from jobs of size <= s
    rho_s = lambda_rate * sympy.integrate(y, (y, 0, s))

    # E_S_s_squared: second moment of service time for jobs of size <= s
    E_S_s_squared = sympy.integrate(y**2, (y, 0, s))

    # 3. Calculate E[W_s], the mean waiting time for a job of size s
    # The correct formula for SRPT has a squared denominator.
    E_W_s_expr = (lambda_rate * E_S_s_squared) / (2 * (1 - rho_s)**2)
    
    # 4. Calculate E[W], the total mean waiting time
    # Integrate E_W_s over the distribution of s (from 0 to 1)
    E_W = sympy.integrate(E_W_s_expr, (s, 0, 1))

    # 5. Calculate x, the optimal mean response time
    x = E_S + E_W

    # 6. Deconstruct x to find the required term
    # x is of the form: rational_part + remaining_part
    # We use sympy to separate the rational part from the rest.
    rational_part = 0
    remaining_term = 0
    
    # x is an Add object in sympy, we can iterate through its arguments
    for term in sympy.Add.make_args(x):
        if term.is_rational:
            rational_part += term
        # Check if the term is a logarithm of a rational number.
        # It must be of the form log(q), not C*log(q).
        elif term.func == sympy.log and term.args[0].is_rational:
             # This case doesn't occur here, but is included for completeness
             rational_part += term
        else:
            remaining_term += term

    # The problem asks for the remaining term after removing rationals and logs of rationals.
    # In our case, x = 7/6 - (2/9)*ln(4).
    # 7/6 is rational and is removed.
    # -(2/9)*ln(4) is not rational and not of the form ln(q), so it remains.
    
    # We need to output each number in the final equation: -2/9, ln(4)
    # Let's get the coefficient and the log term.
    coeff, log_term = remaining_term.as_coeff_Mul()

    # Format for LaTeX
    latex_str = sympy.latex(coeff) + r" \ln(" + sympy.latex(log_term.args[0]) + ")"
    print(latex_str)


solve_queueing_problem()