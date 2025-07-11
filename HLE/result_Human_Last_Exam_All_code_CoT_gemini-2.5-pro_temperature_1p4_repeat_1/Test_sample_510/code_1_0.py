import sympy

def solve_queueing_problem():
    """
    This function calculates the optimal mean response time x and identifies the required term.
    """
    # 1. Define the symbolic variable and system parameters
    s = sympy.Symbol('s')
    lambda_val = sympy.Rational(3, 2)

    # 2. Define terms for the SRPT formula for U(0,1) distribution (f(t)=1)
    
    # rho_s = lambda * integral(t*f(t), (t,0,s))
    # integral(t*1, (t,0,s)) = s**2/2
    rho_s = lambda_val * (s**2 / 2)
    
    # Numerator of the second part of E[T(s)]
    # lambda * integral(t**2*f(t), (t,0,s)) = lambda * s**3/3
    lambda_E_Ss_2 = lambda_val * (s**3 / 3)
    
    # 3. Define the integrand E[T(s)]
    # E[T(s)] consists of two terms. We will integrate them separately.
    term1_integrand = lambda_E_Ss_2 / (2 * (1 - rho_s)**2)
    term2_integrand = s / (1 - rho_s)
    
    # 4. Perform the definite integration from s=0 to s=1
    integral1 = sympy.integrate(term1_integrand, (s, 0, 1))
    integral2 = sympy.integrate(term2_integrand, (s, 0, 1))
    
    x = integral1 + integral2
    
    # 5. Analyze the result and format the final answer
    # The calculated x is 2/3 + 8*ln(2)/9.
    # The additive rational term is 2/3.
    rational_term = sympy.Rational(2, 3)
    
    # The other additive term is 8*ln(2)/9.
    # We check if this term is a logarithm of a rational number.
    # It can be written as ln(2^(8/9)). Since 2^(8/9) is not rational, this term is not removed.
    remaining_term = 8 * sympy.log(2) / 9

    print(f"The optimal mean response time x is the sum of two integrals.")
    print(f"Value of the first integral: {integral1}")
    print(f"Value of the second integral: {integral2}")
    print(f"Total optimal mean response time x = {x}")
    print("-" * 20)
    print(f"The expression for x is composed of an additive rational term '{rational_term}' and another term '{remaining_term}'.")
    print(f"As per the instructions, we remove the additive rational term.")
    print(f"The remaining term is {remaining_term}.")
    
    # For the final equation a/b * ln(c)
    a = 8
    b = 9
    c = 2
    
    print(f"The numbers that form the final expression a/b * ln(c) are: a={a}, b={b}, c={c}.")

if __name__ == '__main__':
    solve_queueing_problem()