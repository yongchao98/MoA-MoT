import math

def solve_integral_approximation():
    """
    Calculates the coefficients for the asymptotic approximation of the integral I(epsilon).
    The approximation is of the form I(epsilon) ~ C1 * epsilon^a + C2 * epsilon^b.
    """
    # Parameters from the integral's denominator polynomial
    # f(x) = p0*x^n + p1*x^k + ...
    p0 = 9.0
    p1 = 5.0
    n = 5.0
    k = 6.0

    # --- Calculate C1, the coefficient of the leading term ---
    # The leading term exponent is a = (1-n)/n
    a = (1 - n) / n
    # The coefficient C1 is given by (1/p0^(1/n)) * integral(1/(1+y^n))
    # The integral is (pi/n) / sin(pi/n)
    integral_c1 = (math.pi / n) / math.sin(math.pi / n)
    c1 = (1 / (p0**(1/n))) * integral_c1

    # --- Calculate C2, the coefficient of the first correction term ---
    # The correction term exponent is b = (k - 2n + 1)/n
    b = (k - 2*n + 1) / n
    # The coefficient C2 involves an integral of the form integral(y^k / (1+y^n)^2)
    # This integral can be simplified to ((k-n+1)/n) * integral(y^(k-n) / (1+y^n))
    # The second integral is (pi/n) / sin(pi*(k-n+1)/n)
    j = k - n
    integral_c2_part2 = (math.pi / n) / math.sin(math.pi * (j + 1) / n)
    integral_c2 = ((j + 1) / n) * integral_c2_part2
    
    # The full C2 coefficient is -p1 / p0^((k+1)/n) * integral_c2
    c2 = - (p1 / (p0**((k+1)/n))) * integral_c2
    
    print("The analytical formula that approximates I(epsilon) for small epsilon is:")
    print(f"I(epsilon) ~= C1 * epsilon^a + C2 * epsilon^b")
    print("\nWhere:")
    print(f"C1 = {c1:.10f}")
    print(f"a = {a}")
    print(f"C2 = {c2:.10f}")
    print(f"b = {b}")
    
    print("\nFinal equation:")
    # Using '+' for the second term and letting the sign of C2 handle it
    print(f"I(epsilon) ~= {c1:.4f} * epsilon**({a:.1f}) + ({c2:.4f}) * epsilon**({b:.1f})")
    print(f"I(epsilon) ~= {c1:.4f} * epsilon**(-4/5) {c2:+.4f} * epsilon**(-3/5)")


solve_integral_approximation()
<<<I(epsilon) ~= 0.6889 * epsilon**(-0.8) - 0.0610 * epsilon**(-0.6)>>>