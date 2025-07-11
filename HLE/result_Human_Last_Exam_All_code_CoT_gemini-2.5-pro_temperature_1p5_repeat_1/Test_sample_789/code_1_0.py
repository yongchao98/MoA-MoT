import math

def get_final_answer():
    """
    This function demonstrates the step-by-step calculation of the rock's mass
    under Titan 5-bit architecture constraints to find the smallest derivable absolute error.
    """

    # --- Step 1: Define initial values as fractions ---
    density_n, density_d = 9, 10
    four_thirds_n, four_thirds_d = 4, 3
    pi_n, pi_d = 22, 7  # 5-bit approximation for pi
    r_n, r_d = 1, 2

    print("Titan-5b Calculation for Mass of a Rock")
    print("Formula: mass = density * (4/3) * pi * r^3")
    print(f"Initial fractions: density={density_n}/{density_d}, pi~={pi_n}/{pi_d}, r={r_n}/{r_d}\n")

    # --- Step 2: Calculate r^3 ---
    r_cubed_n = r_n**3
    r_cubed_d = r_d**3
    print(f"Step 1: Calculate r^3 = ({r_n}/{r_d})^3 = {r_cubed_n}/{r_cubed_d}")
    # This is a valid 5-bit fraction (1 and 8 are <= 31)

    # --- Step 3: Group terms and calculate ---
    # Grouping as [ (9/10)*(22/7) ] * [ (4/3)*(1/8) ]
    print("\nStep 2: Grouping terms as [ (9/10) * (22/7) ] * [ (4/3) * (1/8) ]")

    # Calculate the second group: (4/3) * (1/8)
    term2_n_pre = four_thirds_n * r_cubed_n
    term2_d_pre = four_thirds_d * r_cubed_d
    # Simplify 4/24 to 1/6 by dividing by GCD
    common_divisor_2 = math.gcd(term2_n_pre, term2_d_pre)
    term2_n = term2_n_pre // common_divisor_2
    term2_d = term2_d_pre // common_divisor_2
    print(f"  - Calculating term [ (4/3) * (1/8) ] = {term2_n_pre}/{term2_d_pre}")
    print(f"    Simplifies to {term2_n}/{term2_d} (valid 5-bit fraction)")

    # Calculate the first group: (9/10) * (22/7)
    term1_n_pre = density_n * pi_n
    term1_d_pre = density_d * pi_d
    print(f"\n  - Calculating term [ (9/10) * (22/7) ] = {term1_n_pre}/{term1_d_pre}")
    print(f"    Result is invalid, as numerator {term1_n_pre} > 31.")

    # Approximate the invalid term
    # The actual value is 198/70 ~= 2.828
    # A close 5-bit approximation is 14/5 = 2.8
    approx_term1_n, approx_term1_d = 14, 5
    print(f"    Approximating with a valid 5-bit fraction: {approx_term1_n}/{approx_term1_d}")

    # --- Step 4: Final calculation ---
    print("\nStep 3: Multiplying the results of the two groups")
    final_n_pre = approx_term1_n * term2_n
    final_d_pre = approx_term1_d * term2_d
    print(f"  - mass_approx = ({approx_term1_n}/{approx_term1_d}) * ({term2_n}/{term2_d}) = {final_n_pre}/{final_d_pre}")

    # Simplify final fraction
    common_divisor_final = math.gcd(final_n_pre, final_d_pre)
    final_n = final_n_pre // common_divisor_final
    final_d = final_d_pre // common_divisor_final
    print(f"  - Simplified final result: {final_n}/{final_d}\n")

    # --- Step 5: Error Calculation ---
    print("Step 4: Calculating the absolute error")
    
    # Final derived equation for mass
    print(f"The final equation is: {final_n} / {final_d}")

    # True value using high-precision pi
    true_mass = 0.9 * (4.0/3.0) * math.pi * (0.5**3)
    calculated_mass = final_n / final_d
    absolute_error = abs(true_mass - calculated_mass)

    print(f"  - True mass: ~{true_mass:.5f}")
    print(f"  - Calculated mass: {calculated_mass:.5f}")
    print(f"  - Absolute Error: {absolute_error:.5f}")

    # Rounded to 0.001 as requested
    error_rounded = round(absolute_error, 3)
    print(f"\nSmallest absolute error (e), rounded to 0.001: {error_rounded}")
    return error_rounded

# --- Execute and print the final answer ---
final_error = get_final_answer()
print(f"<<<{final_error}>>>")