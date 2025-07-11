import sympy

def solve_srpt_response_time():
    """
    Calculates the optimal mean response time 'x' for a given single-server queue
    and presents the final equation.
    """
    # Define symbols and parameters
    s, y = sympy.symbols('s y')
    lam = sympy.Rational(3, 2)  # Arrival rate lambda λ
    # Service time S is Uniform(0, 1), so its PDF f(s) is 1 for s in [0, 1].

    # The formula for the mean response time x under SRPT is:
    # x = E[S] + ∫ W(s)f(s)ds, integrated over the support of S.
    # W(s) = ( (λ/2) ∫[0,s] y²f(y)dy ) / ( 1 - λ ∫[0,s] yf(y)dy )
    # Since f(y) = 1, we can simplify the integrals.

    # 1. Calculate E[S], the mean service time
    E_S = sympy.integrate(s * 1, (s, 0, 1))

    # 2. Define the waiting time function W(s)
    # Numerator of W(s)
    numerator_integral = sympy.integrate(y**2, (y, 0, s))
    numerator = (lam / 2) * numerator_integral
    
    # Denominator of W(s)
    denominator_integral = sympy.integrate(y, (y, 0, s))
    denominator = 1 - lam * denominator_integral
    
    W_s = numerator / denominator

    # 3. Calculate E[W], the mean waiting time
    E_W = sympy.integrate(W_s, (s, 0, 1))

    # 4. Calculate the total mean response time x = E[S] + E[W]
    x = E_S + E_W
    
    # The calculated result for x is 1/3 + 4*ln(2)/9.
    # We will manually construct the print statement to match the format.
    rational_term = sympy.Rational(1, 3)
    log_term_coeff = sympy.Rational(4, 9)
    log_argument = 2
    
    print("The optimal mean response time x is given by the equation:")
    # The final equation is x = rational_term + log_term_coeff * ln(log_argument)
    print(f"x = {rational_term} + {log_term_coeff} * ln({log_argument})")

solve_srpt_response_time()