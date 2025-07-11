import math

def solve():
    """
    Solves the Titan computer problem by deriving a calculation for the rock's mass.
    """
    MAX_INT = 31

    # --- Step 1: Define initial formula and values ---
    print("Objective: Calculate Mass = (4/3) * pi * density * radius^3")
    print("-" * 20)

    # Initial values as 5-bit integer fractions
    f_four_thirds = (4, 3)
    f_density = (9, 10)  # 0.9 kg/cm^3
    f_radius = (1, 2)    # 0.5 cm
    
    print(f"Initial fractions:")
    print(f"Constant term = {f_four_thirds[0]}/{f_four_thirds[1]}")
    print(f"Density = {f_density[0]}/{f_density[1]}")
    print(f"Radius = {f_radius[0]}/{f_radius[1]}")

    # --- Step 2: Choose pi and calculate r^3 ---
    # Using pi ~ 22/7 leads to intermediate products like 132/35 or 33/70, which are invalid.
    # We choose pi ~ 25/8 = 3.125, which is less precise but allows for simplification.
    f_pi = (25, 8)
    print(f"Approximating pi ≈ {f_pi[0]}/{f_pi[1]}\n")

    # Calculate r^3
    print("Calculating radius^3:")
    print(f"({f_radius[0]}/{f_radius[1]})^3 = ({f_radius[0]}/{f_radius[1]}) * ({f_radius[0]}/{f_radius[1]}) * ({f_radius[0]}/{f_radius[1]})")
    r_sq_num = f_radius[0] * f_radius[0]
    r_sq_den = f_radius[1] * f_radius[1]
    print(f"Intermediate step: {r_sq_num}/{r_sq_den}") # 1/4
    
    r_cubed_num = r_sq_num * f_radius[0]
    r_cubed_den = r_sq_den * f_radius[1]
    f_r_cubed = (r_cubed_num, r_cubed_den)
    print(f"Final radius^3 = {f_r_cubed[0]}/{f_r_cubed[1]}\n") # 1/8
    
    # --- Step 3: Step-by-step multiplication ---
    # The full calculation is (4/3) * (25/8) * (9/10) * (1/8)
    print("Multiplying terms step-by-step, simplifying to stay within 5-bit limits:")
    
    # Multiply (4/3) * pi
    # (4/3) * (25/8) = 100/24 (invalid). Must simplify first.
    # (4/8) * (25/3) = (1/2) * (25/3)
    num1 = 1 * 25
    den1 = 2 * 3
    # Check intermediate values
    if num1 > MAX_INT or den1 > MAX_INT: raise ValueError("Error")
    term1 = (num1, den1) # 25/6
    print(f"(4/3) * (25/8) simplifies to (1/3) * (25/2) = {term1[0]}/{term1[1]}")
    
    # Multiply by density
    # (25/6) * (9/10) = 225/60 (invalid). Must simplify first.
    # (25/10) * (9/6) = (5/2) * (3/2)
    num2 = 5 * 3
    den2 = 2 * 2
    if num2 > MAX_INT or den2 > MAX_INT: raise ValueError("Error")
    term2 = (num2, den2) # 15/4
    print(f"({term1[0]}/{term1[1]}) * (9/10) simplifies to (5/2) * (3/2) = {term2[0]}/{term2[1]}")
    
    # Multiply by r^3
    # (15/4) * (1/8) = 15/32. The denominator 32 is > 31 and is invalid.
    # We must approximate a term to proceed.
    print(f"Next step: ({term2[0]}/{term2[1]}) * ({f_r_cubed[0]}/{f_r_cubed[1]}) gives 15/32. Denominator is too large.")
    
    # We approximate r^3 = 1/8 (0.125) with a close fraction that allows for simplification.
    # 4/31 (approx 0.129) is a good choice, as the '4' cancels with the denominator of 15/4.
    f_r_cubed_approx = (4, 31)
    print(f"Approximating {f_r_cubed[0]}/{f_r_cubed[1]} with {f_r_cubed_approx[0]}/{f_r_cubed_approx[1]} to enable reduction.")
    
    # Perform final calculation
    # (15/4) * (4/31)
    # The 4s cancel
    final_num = 15
    final_den = 31
    final_fraction = (final_num, final_den)
    print(f"Final calculation: ({term2[0]}/{term2[1]}) * ({f_r_cubed_approx[0]}/{f_r_ cubed_approx[1]}) = {final_fraction[0]}/{final_fraction[1]}\n")
    
    # --- Step 4: Calculate the absolute error ---
    titan_result = final_fraction[0] / final_fraction[1]
    
    # True value using high precision math
    true_volume = (4/3) * math.pi * (0.5**3)
    true_mass = true_volume * 0.9
    
    absolute_error = abs(true_mass - titan_result)
    rounded_error = round(absolute_error, 3)
    
    print("--- Error Analysis ---")
    print(f"The derived mass is {final_fraction[0]}/{final_fraction[1]} ≈ {titan_result:.5f} kg")
    print(f"The true mass is ≈ {true_mass:.5f} kg")
    print(f"The absolute error is |{true_mass:.5f} - {titan_result:.5f}| = {absolute_error:.5f}")
    print(f"The absolute error rounded to 0.001 is {rounded_error}")

    return rounded_error

# Run the calculation and get the final error
final_error_value = solve()

# Print the final answer in the specified format
print(f"<<<{final_error_value}>>>")