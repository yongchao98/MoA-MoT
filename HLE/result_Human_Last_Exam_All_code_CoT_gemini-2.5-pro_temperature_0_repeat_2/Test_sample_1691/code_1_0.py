import math

def solve():
    """
    Calculates the coefficients and exponents for the asymptotic approximation of the integral.
    The integral is I(eps) = integral from 0 to 15 of 1 / (eps + 9x^5 + 5x^6 + 9x^8) dx.
    The approximation has the form I(eps) ~ C1 * eps^e1 + C2 * eps^e2.
    """
    
    # Parameters from the denominator polynomial g(x) = p1*x^a + p2*x^b + p3*x^c
    p1 = 9.0
    a = 5.0
    p2 = 5.0
    b = 6.0
    
    # --- Exponent Calculation ---
    # Exponent of the leading term
    e1 = (1 - a) / a
    # Exponent of the first correction term
    e2 = (b - 2 * a + 1) / a
    
    # --- Coefficient Calculation ---
    pi = math.pi
    
    # Coefficient of the leading term, C1
    # C1 = p1^(-1/a) * integral_0^inf (1 / (1+u^a)) du
    # The integral is pi / (a * sin(pi/a))
    integral_1 = pi / (a * math.sin(pi / a))
    C1 = (p1**(-1/a)) * integral_1
    
    # Coefficient of the second term, C2
    # C2 = -p2 * p1^(-(b+1)/a) * integral_0^inf (u^b / (1+u^a)^2) du
    # The integral can be shown to be ((b-a+1)/a) * integral_0^inf (u^(b-a) / (1+u^a)) du
    # The second integral is pi / (a * sin(pi*(b-a+1)/a))
    integral_2_part1 = (b - a + 1) / a
    integral_2_part2 = pi / (a * math.sin(pi * (b - a + 1) / a))
    integral_2 = integral_2_part1 * integral_2_part2
    C2 = -p2 * (p1**(-(b + 1) / a)) * integral_2
    
    # --- Output the results ---
    print("The analytical approximation for I(epsilon) for small epsilon is of the form:")
    print("I(eps) = C1 * eps^e1 + C2 * eps^e2")
    print("\nCalculated values for the constants:")
    print(f"C1 = {C1}")
    print(f"e1 = {e1}")
    print(f"C2 = {C2}")
    print(f"e2 = {e2}")
    
    print("\nFinal formula with numerical values:")
    # We need to output each number in the final equation.
    # The '+' sign is handled by the sign of C2.
    print(f"I(eps) = {C1} * eps**({e1}) + ({C2}) * eps**({e2})")

solve()