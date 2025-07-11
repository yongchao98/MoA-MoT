import math

def evaluate_integral(zeta_i, zeta_j):
    """
    Evaluates the integral <phi_i | 1/r | phi_j> for 1s Slater orbitals.
    
    The general formula is I = (4 * (zeta_i * zeta_j)**1.5) / (zeta_i + zeta_j)**2.
    For the special case where zeta_i = zeta_j = zeta, the integral simplifies to zeta.
    This function demonstrates the calculation using the general formula.
    """
    
    print(f"Evaluating for zeta_i = {zeta_i} and zeta_j = {zeta_j}\n")
    
    # Numerator calculation
    term1 = zeta_i * zeta_j
    term2 = math.pow(term1, 1.5)
    numerator = 4 * term2
    
    # Denominator calculation
    term3 = zeta_i + zeta_j
    denominator = math.pow(term3, 2)
    
    # Final result
    result = numerator / denominator
    
    # Print the step-by-step evaluation of the formula
    print("The formula is: I = (4 * (zeta_i * zeta_j)**1.5) / (zeta_i + zeta_j)**2\n")
    print("--- Calculation Steps ---")
    print(f"1. zeta_i * zeta_j = {zeta_i} * {zeta_j} = {term1}")
    print(f"2. (zeta_i * zeta_j)**1.5 = ({term1})**1.5 = {term2}")
    print(f"3. Numerator = 4 * {term2} = {numerator}\n")
    
    print(f"4. zeta_i + zeta_j = {zeta_i} + {zeta_j} = {term3}")
    print(f"5. Denominator = (zeta_i + zeta_j)**2 = ({term3})**2 = {denominator}\n")
    
    print(f"Final Result: I = Numerator / Denominator = {numerator} / {denominator} = {result}")
    
    if zeta_i == zeta_j:
        print(f"\nAs expected for this special case, the result is equal to zeta ({zeta_i}).")
        
    return result

# The question asks to evaluate for "a 1s Slater orbital", which implies
# we are finding the expectation value, so phi_i = phi_j.
# This means zeta_i = zeta_j = zeta.
# We will use zeta = 1.24, a common textbook value for an effective orbital exponent.
zeta = 1.24
final_value = evaluate_integral(zeta, zeta)