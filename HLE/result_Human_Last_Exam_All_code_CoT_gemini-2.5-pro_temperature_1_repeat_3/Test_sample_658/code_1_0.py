import math

def print_asymptotic_formula():
    """
    This function constructs and prints the asymptotic formula for d_{B,delta}.
    The formula is A(B, L) = C1 * sqrt(B*L) + C2 * L / log(C3*L/B + C4),
    where L = log(1/delta).

    The asymptotic analysis reveals the following constants:
    """
    
    # The coefficient for the first term, sqrt(B*L)
    C1 = 1
    
    # The coefficient for L in the numerator of the second term
    C2 = 1
    
    # The coefficient for the ratio L/B inside the logarithm
    C3 = 1
    
    # A constant added to the argument of the logarithm.
    # It must be greater than 1 to ensure the argument is always > 1, since L/B can approach 0.
    # The value of this constant does not affect the Theta-asymptotics. We choose 2.
    C4 = 2

    # We construct the formula string part-by-part to handle coefficients equal to 1.
    
    term1 = f"sqrt(B*L)"
    if C1 != 1:
        term1 = f"{C1}*{term1}"
        
    term2_numerator = "L"
    if C2 != 1:
        term2_numerator = f"{C2}*L"
        
    term2_denominator_arg = "L/B"
    if C3 != 1:
        term2_denominator_arg = f"{C3}*L/B"
        
    term2_denominator = f"log({term2_denominator_arg} + {C4})"
    
    term2 = f"{term2_numerator} / {term2_denominator}"
    
    final_formula = f"{term1} + {term2}"
    
    print("The asymptotic value A(B, L) for d_{B,delta}, where L = log(delta^-1), is given by:")
    print(final_formula)
    print("\nThe numbers (constants) in this formula are:")
    print(f"C1 = {C1} (for the first term)")
    print(f"C2 = {C2} (for the numerator of the second term)")
    print(f"C3 = {C3} (for the L/B ratio in the log)")
    print(f"C4 = {C4} (the constant added in the log)")

print_asymptotic_formula()