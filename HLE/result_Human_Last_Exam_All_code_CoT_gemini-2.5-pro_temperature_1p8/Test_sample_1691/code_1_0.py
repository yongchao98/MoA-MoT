import math

def approximate_integral_formula():
    """
    This function derives and prints the analytical formula that approximates
    I(epsilon) for small epsilon.
    """
    
    # The integral is of the form integral(1/(epsilon + f(x))) dx
    # where f(x) = 9.0*x**5.0 + 5.0*x**6.0 + 9.0*x**8.0.
    # For small x, the behavior of f(x) is dominated by the term with the
    # lowest power, which is a*x**p.
    a = 9.0
    p = 5.0
    
    # The approximation for small epsilon is of the form I(epsilon) ≈ C * epsilon**k.
    
    # The exponent k is derived from a scaling analysis and is equal to -(p-1)/p.
    k_numerator = int(-(p - 1))
    k_denominator = int(p)
    k_val = k_numerator / k_denominator
    
    # The coefficient C is given by the formula: pi / (p * a**(1/p) * sin(pi/p)).
    C_val = math.pi / (p * (a**(1.0/p)) * math.sin(math.pi/p))
    
    # Print the detailed analytical formula.
    print("The analytical formula that approximates I(epsilon) for the small epsilon regime is:")
    # We display each number that constitutes the final equation as requested.
    print(f"I(epsilon) ≈ (π / ({p} * ({a})**(1/{p}) * sin(π/{p}))) * epsilon**({k_numerator}/{k_denominator})")
    
    # Print the numerically evaluated version of the formula for practical use.
    print("\nThis evaluates to the approximate formula:")
    print(f"I(epsilon) ≈ {C_val:.5f} * epsilon**({k_val})")

if __name__ == '__main__':
    approximate_integral_formula()
