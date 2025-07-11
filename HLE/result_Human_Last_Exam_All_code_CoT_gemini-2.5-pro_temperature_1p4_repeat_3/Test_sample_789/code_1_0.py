import math

def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    while b:
        a, b = b, a % b
    return a

def multiply_fractions(f1, f2):
    """
    Multiplies two fractions (n1/d1) * (n2/d2) under Titan's 5-bit constraint.
    Simplifies before multiplying to prevent overflow.
    """
    n1, d1 = f1
    n2, d2 = f2
    
    # Check for invalid inputs (though not expected in this script's flow)
    if any(x > 31 for x in [n1, d1, n2, d2]):
        raise ValueError("Input fraction contains number > 31")

    # Simplify n1 with d2
    common1 = gcd(n1, d2)
    n1_s = n1 // common1
    d2_s = d2 // common1
    
    # Simplify n2 with d1
    common2 = gcd(n2, d1)
    n2_s = n2 // common2
    d1_s = d1 // common2
    
    # Perform multiplication on simplified parts
    num = n1_s * n2_s
    den = d1_s * d2_s
    
    # Check for overflow in the result
    if num > 31 or den > 31:
        raise ValueError(f"Multiplication resulted in overflow: {num}/{den}")
        
    return (num, den)

def main():
    """
    Derives the mass calculation and computes the smallest absolute error.
    """
    # Problem values
    r = 0.5
    rho = 0.9

    # --- True Mass Calculation (for comparison) ---
    true_mass = (4/3) * math.pi * (r**3) * rho

    # --- Titan 5-bit Fractional Calculation ---
    print("--- Titan 5-bit Calculation ---")
    print("Goal: Calculate mass m = (4/3) * pi * r^3 * rho with minimal error.")
    
    # Step 1: Define inputs as fractions
    f_four_thirds = (4, 3)
    f_r_cubed = (1, 8) # r = 1/2, so r^3 = 1/8
    
    # Best 5-bit approximation for pi
    f_pi = (22, 7)
    
    # The exact fraction for rho = 0.9 is 9/10.
    # A direct calculation with 9/10 leads to overflow (3*11=33).
    # We approximate rho = 0.9 with 10/11 = 0.909... to enable simplification.
    f_rho_approx = (10, 11)
    
    print("\nInitial equation with chosen fractions:")
    print(f"m = ({f_four_thirds[0]}/{f_four_thirds[1]}) * ({f_pi[0]}/{f_pi[1]}) * ({f_r_cubed[0]}/{f_r_cubed[1]}) * ({f_rho_approx[0]}/{f_rho_approx[1]})")

    # Step 2: Perform the calculation, showing intermediate steps
    print("\nCalculation steps (simplifying to avoid overflow):")
    
    # Combine pi and rho terms: (22/7) * (10/11)
    # This simplifies to (2/7) * (10/1) = 20/7
    term1 = multiply_fractions(f_pi, f_rho_approx)
    print(f"1. Combine pi and density: ({f_pi[0]}/{f_pi[1]}) * ({f_rho_approx[0]}/{f_rho_approx[1]}) = {term1[0]}/{term1[1]}")
    
    # Combine the constant and radius terms: (4/3) * (1/8)
    # This simplifies to (1/3) * (1/2) = 1/6
    term2 = multiply_fractions(f_four_thirds, f_r_cubed)
    print(f"2. Combine constant and radius: ({f_four_thirds[0]}/{f_four_thirds[1]}) * ({f_r_cubed[0]}/{f_r_cubed[1]}) = {term2[0]}/{term2[1]}")

    # Combine the intermediate results: (1/6) * (20/7)
    # This simplifies to (1/3) * (10/7) = 10/21
    final_mass_frac = multiply_fractions(term2, term1)
    print(f"3. Combine results: ({term2[0]}/{term2[1]}) * ({term1[0]}/{term1[1]}) = {final_mass_frac[0]}/{final_mass_frac[1]}")
    
    # Step 3: Final result and error calculation
    calculated_mass = final_mass_frac[0] / final_mass_frac[1]
    
    print("\n--- Results ---")
    print(f"Final calculated mass: {final_mass_frac[0]}/{final_mass_frac[1]} ≈ {calculated_mass:.5f} kg")
    print(f"True mass:             ≈ {true_mass:.5f} kg")

    abs_error = abs(true_mass - calculated_mass)
    rounded_error = round(abs_error, 3)

    print(f"Smallest absolute error (e): {abs_error:.5f}")
    print(f"Error rounded to 0.001 is: {rounded_error}")
    
if __name__ == "__main__":
    main()
