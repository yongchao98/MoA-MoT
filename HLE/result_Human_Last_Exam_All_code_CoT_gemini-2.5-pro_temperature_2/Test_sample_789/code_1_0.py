import math

def calculate_mass_on_titan():
    """
    Performs the mass calculation following Titan 5-bit architecture rules.
    """
    
    # Helper function to find the greatest common divisor for simplification
    def gcd(a, b):
        while b:
            a, b = b, a % b
        return a

    print("Step 1: Define initial values as 5-bit fractions.")
    # Density: 0.9 kg/cm^3 -> 9/10
    density = (9, 10)
    # Volume constant from formula
    vol_const = (4, 3)
    # Pi approximation (best 5-bit fraction)
    pi = (22, 7)
    # Radius: 0.5 cm -> 1/2
    radius = (1, 2)
    # Radius cubed: (1/2)^3 -> 1/8
    r_cubed = (radius[0]**3, radius[1]**3)

    print(f"  Density (d)       = {density[0]}/{density[1]}")
    print(f"  Volume Constant   = {vol_const[0]}/{vol_const[1]}")
    print(f"  Pi Approximation  = {pi[0]}/{pi[1]}")
    print(f"  Radius Cubed (r^3) = {r_cubed[0]}/{r_cubed[1]}\n")

    print("Step 2: Combine terms step-by-step, simplifying to stay within 5-bits.")
    # Full expression: Mass = (9/10) * (4/3) * (22/7) * (1/8)
    # We combine (9/10), (4/3), and (1/8) first.
    # Let's multiply them together: (9*4*1) / (10*3*8) = 36 / 240
    combined_num = density[0] * vol_const[0] * r_cubed[0]
    combined_den = density[1] * vol_const[1] * r_cubed[1]
    
    # Simplify the combined fraction
    common_divisor = gcd(combined_num, combined_den)
    constants_frac = (combined_num // common_divisor, combined_den // common_divisor)
    print(f"  Combining non-pi terms: (9/10 * 4/3 * 1/8) = 36/240 simplifies to {constants_frac[0]}/{constants_frac[1]}.")
    
    print("\nStep 3: Multiply by the Pi approximation.")
    # Expression is now M = (3/20) * (22/7)
    print(f"  The expression becomes: M = ({constants_frac[0]}/{constants_frac[1]}) * ({pi[0]}/{pi[1]})")

    # This multiplication results in (3 * 22) / (20 * 7) = 66 / 140.
    # Simplified, this is 33/70.
    invalid_num = 33
    invalid_den = 70
    print(f"  Multiplying gives {invalid_num}/{invalid_den}. The numerator and denominator are larger than 31 and thus invalid.\n")

    print("Step 4: Replace the invalid fraction with its closest 5-bit approximation.")
    target_value = invalid_num / invalid_den
    print(f"  The target value is {invalid_num}/{invalid_den} ≈ {target_value:.5f}")
    
    best_frac = (0, 1)
    min_diff = float('inf')

    # Search all possible 5-bit fractions for the best fit
    for a in range(1, 32):
        for b in range(1, 32):
            val = a / b
            diff = abs(val - target_value)
            if diff < min_diff:
                min_diff = diff
                best_frac = (a, b)
    
    final_mass_frac = best_frac
    print(f"  The best 5-bit fractional replacement is {final_mass_frac[0]}/{final_mass_frac[1]}.\n")

    print("Step 5: State final derived equation and calculate the error.")
    print("The full calculation is:")
    final_equation_str = (f"  {density[0]}/{density[1]} * {vol_const[0]}/{vol_const[1]} * "
                          f"{pi[0]}/{pi[1]} * {r_cubed[0]}/{r_cubed[1]} = {final_mass_frac[0]}/{final_mass_frac[1]}")
    print(final_equation_str)

    # Calculate final numerical values
    final_mass_val = final_mass_frac[0] / final_mass_frac[1]
    true_mass = 0.9 * (4/3) * math.pi * (0.5**3) # approx 0.15 * pi
    absolute_error = abs(final_mass_val - true_mass)
    
    print(f"\n  Derived Mass   = {final_mass_val:.7f} kg")
    print(f"  True Mass      = {true_mass:.7f} kg")
    print(f"  Absolute Error = {absolute_error:.7f}\n")
    
    # Round the absolute error to 0.001 as required by the problem
    error_rounded = round(absolute_error * 1000) / 1000.0
    
    print("Final answer (e), which is the smallest absolute error, rounded to 0.001.")
    
    return error_rounded

# --- Execute the calculation ---
if __name__ == "__main__":
    final_error = calculate_mass_on_titan()
    print(f"<<<{final_error}>>>")