import sympy

def solve_integral():
    """
    This function calculates the value of the complex integral provided by the user.
    The steps are outlined in the thought process.
    """
    pi = sympy.pi
    p = sympy.symbols('p')

    # Step 1: Denominator simplification.
    # The maximization problem is max(4 * sum((exp(m_i)-1)**2)) for sum(m_i**2) <= p**2.
    # The maximum value is 4 * (exp(p) - 1)**2 for p > 0.
    # The denominator is sqrt(4 * (exp(p) - 1)**2) = 2 * (exp(p) - 1).

    # Step 2: Simplify the integrand.
    # Numerator: 2*p - exp(-p/4) + 2*p**7 + 2*p*exp(-p) + exp(p/4)
    # Denominator: 2*(exp(p) - 1)
    # Simplified integrand: (p + p**7 + p*exp(-p) + sinh(p/4)) / (exp(p) - 1)
    # We split this into four integrals I1, I2, I3, I4.

    # Step 3 & 4: Calculate each integral.

    # I1 = integral from 0 to inf of p / (exp(p) - 1) dp
    # This is Gamma(2) * zeta(2)
    s1 = 2
    I1 = sympy.gamma(s1) * sympy.zeta(s1)
    # I1 = 1 * pi**2 / 6

    # I2 = integral from 0 to inf of p**7 / (exp(p) - 1) dp
    # This is Gamma(8) * zeta(8)
    s2 = 8
    I2 = sympy.gamma(s2) * sympy.zeta(s2)
    # Gamma(8) = 7! = 5040. zeta(8) = pi**8/9450.
    # I2 = 5040 * pi**8 / 9450 = (504/945) * pi**8 = (8/15) * pi**8

    # I3 = integral from 0 to inf of p*exp(-p) / (exp(p) - 1) dp
    # We rewrite the integrand part as p * (1/(exp(p)-1) - exp(-p))
    # I3 = Integral(p/(exp(p)-1) dp) - Integral(p*exp(-p) dp)
    # The first part is I1. The second part is Gamma(2) = 1.
    I3 = I1 - sympy.gamma(2)

    # I4 = integral from 0 to inf of sinh(p/4) / (exp(p) - 1) dp
    # This evaluates to 2 - pi/2 via a series sum.
    # The derivation leads to the sum (1/4) * Sum_{k=1 to inf} 1/(k**2 - (1/4)**2)
    # Using the formula Sum_{k=1 to inf} 1/(k**2 - a**2) = 1/(2a**2) - (pi*cot(pi*a))/(2a)
    # For a=1/4, the sum is 8 - 2*pi.
    # I4 = (1/4) * (8 - 2*pi) = 2 - pi/2.
    a = sympy.Rational(1, 4)
    # Not using the formula directly, just the result of the calculation.
    I4 = 2 - pi / 2

    # Step 5: Sum the results
    total_value = I1 + I2 + I3 + I4
    
    # Let's collect the terms for the final printout.
    term_pi8 = I2
    term_pi2 = I1 + I3.coeff(pi**2) # I1 + I1 from I3
    term_pi = I4.coeff(pi) * pi
    const_term = I3.as_coeff_add()[0] + I4.as_coeff_add()[0] # -1 + 2

    print("The final expression is a sum of four calculated integrals:")
    print(f"Integral of p^7 term: {term_pi8}")
    print(f"Integral of p term: {I1}")
    print(f"Integral of p*exp(-p) term: {I3}")
    print(f"Integral of sinh(p/4) term: {I4}")
    
    print("\nCombining these terms gives the total value:")
    # Using sympy.collect to group terms by powers of pi
    final_expr = sympy.collect(total_value, [pi**8, pi**2, pi])
    
    # Extract coefficients for printing the final equation
    c_pi8 = final_expr.coeff(pi**8)
    c_pi2 = final_expr.coeff(pi**2)
    c_pi = final_expr.coeff(pi)
    const = final_expr.subs(pi, 0)
    
    # Print the equation with each number explicitly shown
    print(f"Final Equation: ({c_pi8})*pi^8 + ({c_pi2})*pi^2 + ({c_pi})*pi + ({const})")

solve_integral()