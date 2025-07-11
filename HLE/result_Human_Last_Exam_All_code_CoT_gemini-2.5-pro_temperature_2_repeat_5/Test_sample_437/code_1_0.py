import math

def evaluate_1s_integral(zeta_i, zeta_j):
    """
    Calculates the integral <phi_i | 1/r | phi_j> for 1s Slater orbitals.

    This function first prints a detailed explanation of the derivation and then
    calculates the integral value based on the provided orbital exponents.

    Args:
        zeta_i (float): Orbital exponent for the first 1s orbital (phi_i).
        zeta_j (float): Orbital exponent for the second 1s orbital (phi_j).

    Returns:
        float: The value of the integral.
    """
    
    # --- Part 1: Explanation and Derivation ---
    print("Evaluating the integral: <φ_i| 1/r |φ_j> for 1s Slater-Type Orbitals (STOs).")
    print("\nThe normalized 1s STO is given by the formula: φ(r) = (ζ³/π)^(1/2) * exp(-ζr)")
    print("where ζ is the orbital exponent.")
    print("\nThe integral in spherical coordinates is I = ∫∫∫ φ_i*(r) * (1/r) * φ_j(r) * r²sin(θ)drdθdφ")
    print("After integrating the angular part (which gives 4π) and simplifying, we get a radial integral:")
    print("I = 4 * (ζ_i*ζ_j)^(3/2) * ∫[0 to ∞] r * exp(-(ζ_i + ζ_j)r) dr")
    print("\nSolving the radial integral (∫x*exp(-ax)dx = 1/a²) gives the final general formula:")
    print("I = 4 * (ζ_i * ζ_j)^(3/2) / (ζ_i + ζ_j)²")
    print("-" * 60)
    
    # --- Part 2: Numerical Calculation ---
    
    zi_zj_prod = zeta_i * zeta_j
    zi_plus_zj = zeta_i + zeta_j
    
    # Using try-except for potential math domain errors with negative inputs
    try:
        numerator = 4 * math.pow(zi_zj_prod, 1.5)
    except ValueError:
        print("Error: The product of exponents must be non-negative.")
        return None
        
    denominator = math.pow(zi_plus_zj, 2)
    
    if denominator == 0:
        print("Error: The sum of exponents cannot be zero.")
        return None
        
    result = numerator / denominator
    
    print(f"Calculation for exponents ζ_i = {zeta_i} and ζ_j = {zeta_j}:")
    
    # As requested, output each number in the final equation.
    print(f"\nNumerator = 4 * (ζ_i * ζ_j)^(3/2)")
    print(f"          = 4 * ({zeta_i} * {zeta_j})^(1.5)")
    print(f"          = 4 * ({zi_zj_prod})^(1.5)")
    print(f"          = 4 * {math.pow(zi_zj_prod, 1.5)}")
    print(f"          = {numerator}")
    
    print(f"\nDenominator = (ζ_i + ζ_j)²")
    print(f"            = ({zeta_i} + {zeta_j})²")
    print(f"            = ({zi_plus_zj})²")
    print(f"            = {denominator}")
    
    print(f"\nResult = Numerator / Denominator")
    print(f"       = {numerator} / {denominator}")
    print(f"       = {result}")
    
    print("-" * 60)
    return result

# --- Main Execution ---

# Case 1: Identical orbitals (the most common case, representing an expectation value)
# We use the standard exponent for a Hydrogen atom's 1s orbital, ζ = 1.0 (in atomic units).
print("--- Case 1: Identical Orbitals (φ_i = φ_j) ---")
print("This corresponds to calculating the expectation value <φ|1/r|φ>.\n")
zeta_val = 1.0
result1 = evaluate_1s_integral(zeta_val, zeta_val)
print(f"Final answer for ζ_i = ζ_j = {zeta_val}: {result1}")
print("\nNote: When ζ_i = ζ_j = ζ, the formula simplifies from 4*(ζ²)^(3/2)/(2ζ)² to ζ.")
print(f"The result is simply the value of ζ, which is {zeta_val}, matching the calculation.")

# Case 2: An example with different orbitals to show the general formula's use.
# print("\n\n" + "="*60 + "\n")
# print("--- Case 2: Different Orbitals (φ_i ≠ φ_j) ---")
# zeta_1 = 1.0
# zeta_2 = 1.5
# result2 = evaluate_1s_integral(zeta_1, zeta_2)
# print(f"Final answer for ζ_i = {zeta_1}, ζ_j = {zeta_2}: {result2}")
