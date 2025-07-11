import math

def titan_calculation_and_error():
    """
    This function demonstrates the step-by-step calculation of the rock's mass
    under the Titan 5-bit architecture constraints and calculates the final error.
    """
    # --- Initial Values and Setup ---
    # Representing inputs as fractions of 5-bit integers (<= 31)
    rho_n, rho_d = 9, 10  # density = 0.9 = 9/10
    r_n, r_d = 1, 2      # radius = 0.5 = 1/2
    pi_n, pi_d = 22, 7   # pi approximation = 22/7
    c_n, c_d = 4, 3      # constant = 4/3

    print("Titan 5-bit Architecture Calculation")
    print("------------------------------------")
    print("Objective: Calculate Mass = density * (4/3) * pi * radius^3")
    print(f"Initial fractions: density = {rho_n}/{rho_d}, pi ~ {pi_n}/{pi_d}, radius = {r_n}/{r_d}\n")

    # --- Step 1: Calculate r^3 ---
    r3_n = r_n**3
    r3_d = r_d**3
    print(f"Step 1: Calculate radius^3 = ({r_n}/{r_d})^3 = {r3_n}/{r3_d}")

    full_equation = f"Mass = ({rho_n}/{rho_d}) * ({c_n}/{c_d}) * ({pi_n}/{pi_d}) * ({r3_n}/{r3_d})"
    print(f"The full equation is: {full_equation}\n")

    # --- Step 2: Reorder and calculate in parts ---
    print("Step 2: Reorder to manage intermediate values: Mass = [({rho_n}/{rho_d})*({pi_n}/{pi_d})] * [({c_n}/{c_d})*({r3_n}/{r3_d})]\n")

    # --- Step 3: Calculate the second part ---
    part2_n = c_n * r3_n
    part2_d = c_d * r3_d
    # Simplify 4/24 to 1/6
    final_part2_n, final_part2_d = 1, 6
    print(f"Step 3: Calculate second part: ({c_n}/{c_d}) * ({r3_n}/{r3_d}) = {part2_n}/{part2_d}")
    print(f"         Simplifying {part2_n}/{part2_d} gives {final_part2_n}/{final_part2_d}. This is a valid intermediate result.\n")

    # --- Step 4: Handle the first part with approximation ---
    part1_num = rho_n * pi_n
    part1_den = rho_d * pi_d
    print(f"Step 4: Calculate first part: ({rho_n}/{rho_d}) * ({pi_n}/{pi_d}) = {part1_num}/{part1_den}")
    print(f"         Problem: {part1_num} and {part1_den} are > 31, so this is an invalid operation.")
    
    # Justification for approximation
    approx_val = part1_num / part1_den
    approx_n, approx_d = 14, 5
    print(f"         We must approximate the result ({approx_val:.4f}).")
    print(f"         To ensure the next multiplication ({approx_n}/{approx_d} * {final_part2_n}/{final_part2_d}) is valid, the denominator must be <= 5.")
    print(f"         The best approximation with this constraint is {approx_n}/{approx_d} = {approx_n/approx_d}.")
    print(f"         So, we replace the result of the first part with {approx_n}/{approx_d}.\n")

    # --- Step 5: Final Calculation ---
    final_mass_n = approx_n * final_part2_n
    final_mass_d = approx_d * final_part2_d
    
    # Simplify 14/30 to 7/15
    simplified_mass_n, simplified_mass_d = 7, 15
    
    print("Step 5: Perform the final multiplication with the approximated value.")
    final_equation = f"Mass ~ ({approx_n} / {approx_d}) * ({final_part2_n} / {final_part2_d}) = {final_mass_n} / {final_mass_d}"
    print(f"         {final_equation}")
    print(f"         The final simplified result is {simplified_mass_n}/{simplified_mass_d}.\n")

    # --- Step 6: Error Calculation ---
    calculated_mass = simplified_mass_n / simplified_mass_d
    actual_mass = (9/10) * (4/3) * math.pi * (1/2)**3
    abs_error = abs(calculated_mass - actual_mass)
    
    print("--- Error Calculation ---")
    print(f"Final calculated mass as a decimal: {calculated_mass:.6f}")
    print(f"Actual mass (using high-precision pi): {actual_mass:.6f}")
    print(f"Absolute error: |{calculated_mass:.6f} - {actual_mass:.6f}| = {abs_error:.6f}")
    
    # --- Step 7: Final Answer ---
    rounded_error = round(abs_error, 3)
    print(f"The smallest absolute error, rounded to 0.001, is: {rounded_error}")

titan_calculation_and_error()
print("<<<0.005>>>")