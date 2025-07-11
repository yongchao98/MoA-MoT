import math

def solve():
    """
    This function calculates the coefficients for the asymptotic expansion of the integral I(epsilon).
    """
    
    # Coefficients and powers from the polynomial in the denominator
    # P(x) = p1*x**n1 + p2*x**n2 + p3*x**n3
    p1 = 9.0
    n1 = 5.0
    p2 = 5.0
    n2 = 6.0
    p3 = 9.0
    n3 = 8.0

    # The powers of epsilon in the expansion are determined by the powers of x in P(x).
    # The first term in the expansion has power -(n1-1)/n1
    # The second term has power -(n1-2)/n1 (since n2 = n1 + 1)
    power1 = -(n1 - 1.0) / n1
    power2 = -(n1 - 2.0) / n1

    # Calculate the coefficient C1 for the leading term
    # C1 = (pi / n1) / (p1**(1/n1) * sin(pi/n1))
    C1 = (math.pi / n1) / (p1**(1.0/n1) * math.sin(math.pi/n1))
    
    # Calculate the coefficient C2 for the second term
    # The formula is derived from the change-of-variable method explained in the plan.
    # C2 = -2 * (p2/p1) * (1/p1**(2/n1)) * (pi / n1) / sin(2*pi/n1)
    # This simplifies to:
    # C2 = -2 * pi * p2 / (n1 * p1**(1 + 2/n1) * sin(2*pi/n1))
    # From the detailed derivation: C2 = 2*d2 * (pi / (n1*sin(2*pi/n1)))
    # where d2 = -c2/c1**3, c1 = p1**(1/n1), c2 = p2/(n1*p1**(1-1/n1))
    # This leads to C2 = -2 * (p2 * pi) / (n1 * p1**((n1+2)/n1) * math.sin(2*math.pi/n1))

    # Based on the manual derivation: 
    # C2_coeff = 2*d2 * (pi / (n1*sin(2*pi/n1)))
    # d2 = -9**(-7/5)
    C2 = -2 * p1**(- (n1+2)/n1) * (p2/n1) * (math.pi / math.sin(2*math.pi/n1))
    
    # A simpler way to write the coefficient C2 is based on the derivation:
    # C2 = -2 * p1**(-7.0/5.0) * (math.pi/(n1 * math.sin(2*math.pi/n1)))
    # The coefficient of x**(n1+1) is p2.
    # The general formula for the coefficient involves 2*d2, where d2 = -c2/c1**3
    # c1 = p1**(1/n1), c2 = p2 / (n1 * p1**((n1-1)/n1))
    # d2 = -(p2 / (n1 * p1**((n1-1)/n1))) / p1**(3/n1) = -p2 / (n1 * p1**((n1+2)/n1))
    # Term is 2 * d2 * Integral(...)
    # Integral is (pi / (n1*sin(2pi/n1))) * eps**((2-n1)/n1)
    # So C2 = 2 * (-p2 / (n1 * p1**((n1+2)/n1))) * (pi / (n1 * math.sin(2*math.pi/n1)))
    # Let's recalculate and simplify
    term_in_d2 = p2/n1
    d2_denom_p1_power = (n1-1)/n1 + 3/n1
    
    # Let's stick to the manually verified coefficients from the plan
    # C1 = (pi / (5 * 9**(1/5) * sin(pi/5)))
    # C2 = -2 * 9**(-7/5) * (pi / (5 * sin(2*pi/5)))
    
    C1_val = (math.pi / (n1 * p1**(1.0/n1) * math.sin(math.pi/n1)))
    C2_val = -2 * p1**(- (n1+2.0)/n1) * p2 * (math.pi / n1) / math.sin(2.0*math.pi/n1)

    # Let's output the values and the formula
    print(f"The integral can be approximated by the formula: I(eps) = C1 * eps**(-4/5) + C2 * eps**(-3/5)")
    print(f"where the coefficients are:")
    print(f"C1 = {C1_val:.10f}")
    print(f"C2 = {C2_val:.10f}")
    print("\nThe final formula with the numerical coefficients is:")
    final_formula = f"I(epsilon) \u2248 {C1_val:.4f} * epsilon**(-{n1-1:.1f}/{n1:.1f}) + ({C2_val:.4f}) * epsilon**(-{n1-2:.1f}/{n1:.1f})"
    print(final_formula)

solve()