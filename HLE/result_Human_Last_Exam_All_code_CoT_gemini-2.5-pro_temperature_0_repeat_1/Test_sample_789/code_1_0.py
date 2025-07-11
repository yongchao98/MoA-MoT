import math

def solve_titan_mass_calculation():
    """
    Performs the mass calculation following Titan 5-bit architecture rules
    and prints the step-by-step derivation.
    """
    # --- Step 1: Define initial values as 5-bit fractions ---
    rho_n, rho_d = 9, 10
    r_n, r_d = 1, 2
    pi_n, pi_d = 22, 7
    ft_n, ft_d = 4, 3

    print("Derivation of the mass calculation:")
    print("Formula: Mass = density * (4/3) * pi * radius^3")
    print(f"Initial values: density = {rho_n}/{rho_d}, 4/3 = {ft_n}/{ft_d}, pi ≈ {pi_n}/{pi_d}, radius = {r_n}/{r_d}")
    print("-" * 40)

    # --- Step 2: Calculate r^3 ---
    r3_n, r3_d = r_n**3, r_d**3
    print(f"Step 1: Calculate radius^3 = ({r_n}/{r_d})^3 = {r3_n}/{r3_d}")
    print(f"Full equation: Mass = ({rho_n}/{rho_d}) * ({ft_n}/{ft_d}) * ({pi_n}/{pi_d}) * ({r3_n}/{r3_d})")
    print("-" * 40)

    # --- Step 3: Combine terms sequentially, handling overflows ---
    # Combine (9/10) * (4/3) -> 36/30 (overflow) -> simplify to 6/5
    term1_n, term1_d = 6, 5
    print(f"Step 2: Combine ({rho_n}/{rho_d}) * ({ft_n}/{ft_d}).")
    print(f"  Result (9*4)/(10*3) = 36/30 overflows. Simplify by cancellation: (9/3)*(4/10) = 3 * 2/5 = {term1_n}/{term1_d}.")
    
    # Combine (22/7) * (1/8) -> 22/56 (overflow) -> simplify to 11/28
    term2_n, term2_d = 11, 28
    print(f"Step 3: Combine ({pi_n}/{pi_d}) * ({r3_n}/{r3_d}).")
    print(f"  Result (22*1)/(7*8) = 22/56 overflows. Simplify by dividing by 2: {term2_n}/{term2_d}.")
    print(f"Equation is now: Mass = ({term1_n}/{term1_d}) * ({term2_n}/{term2_d})")
    print("-" * 40)

    # --- Step 4: Handle the final multiplication which overflows ---
    print(f"Step 4: Calculate ({term1_n}/{term1_d}) * ({term2_n}/{term2_d}).")
    print(f"  Direct multiplication (6*11)/(5*28) = 66/140 overflows.")
    print(f"  Rearranging to (6/28)*(11/5) and simplifying gives (3/14)*(11/5).")
    print(f"  This multiplication, (3*11)/(14*5) = 33/70, also overflows.")
    print("-" * 40)

    # --- Step 5: Approximate to resolve the overflow ---
    approx_n, approx_d = 13, 6
    print(f"Step 5: Approximate to proceed. We replace 11/5 (value 2.2) with {approx_n}/{approx_d} (value ~2.167).")
    print(f"  The calculation becomes: (3/14) * ({approx_n}/{approx_d})")
    
    # --- Step 6: Final calculation with simplification ---
    final_n, final_d = 13, 28
    print(f"  Multiplying (3*13)/(14*6) = 39/84 overflows. We simplify first:")
    print(f"  (3/6) * (13/14) = (1/2) * (13/14) = {final_n}/{final_d}")
    print("-" * 40)

    # --- Step 7: Calculate and report the final error ---
    true_value = 0.9 * (4.0/3.0) * math.pi * (0.5**3)
    approx_value = final_n / final_d
    absolute_error = abs(true_value - approx_value)
    final_answer = round(absolute_error, 3)

    print("Final Result and Error:")
    print(f"The final calculated mass is {final_n}/{final_d} ≈ {approx_value:.5f} kg.")
    print(f"The true mass is ≈ {true_value:.5f} kg.")
    print(f"The absolute error is |{true_value:.5f} - {approx_value:.5f}| = {absolute_error:.5f}")
    print(f"The smallest absolute error, rounded to 0.001, is {final_answer}.")

    print(f"<<<{final_answer}>>>")

solve_titan_mass_calculation()