import math

def solve_titan_mass_problem():
    """
    This script calculates the mass of a rock based on Titan's 5-bit fractional architecture.
    It follows the rules of the architecture, showing simplifications and approximations.
    """
    # Step 1: Define initial values as 5-bit fractions
    # density, ρ = 0.9 kg/cm³ -> 9/10
    rho_num, rho_den = 9, 10
    # Constant 4/3
    four_thirds_num, four_thirds_den = 4, 3
    # Pi approximation π ≈ 22/7
    pi_num, pi_den = 22, 7
    # Radius, r = 0.5 cm -> 1/2
    r_num, r_den = 1, 2

    print("--- Titan 5-Bit Calculation for Mass of a Rock ---")
    print("\nInitial Formula:")
    print(f"m = density * (4/3) * pi * r^3")
    print(f"m = ({rho_num}/{rho_den}) * ({four_thirds_num}/{four_thirds_den}) * ({pi_num}/{pi_den}) * ({r_num}/{r_den})^3")
    print("-" * 50)

    # --- Calculation Steps ---
    print("Calculation Steps:")

    # Step A: Calculate r³
    r_cubed_num = r_num**3
    r_cubed_den = r_den**3
    print(f"1. Calculate r^3: ({r_num}/{r_den})^3 = {r_cubed_num}/{r_cubed_den}")

    # Step B: Multiply (pi) * (r³)
    # Direct multiplication would be (22*1)/(7*8) = 22/56. Denominator 56 > 31 (Overflow).
    # We must simplify this result to proceed.
    res1_num_unsimplified = pi_num * r_cubed_num
    res1_den_unsimplified = pi_den * r_cubed_den
    # Simplify 22/56 by dividing by gcd(22, 56)=2. Result is 11/28.
    res1_num = 11
    res1_den = 28
    print(f"2. Multiply by pi: ({pi_num}/{pi_den}) * ({r_cubed_num}/{r_cubed_den}) = {res1_num_unsimplified}/{res1_den_unsimplified} -> OVERFLOW")
    print(f"   Simplifying gives: {res1_num}/{res1_den} (Valid)")

    # Step C: Multiply (4/3) * result from Step B
    # Direct multiplication would be (4*11)/(3*28) = 44/84. Both > 31 (Overflow).
    # We simplify before multiplying by cancelling the common factor of 4.
    # (4/3) * (11/28) -> (1/3) * (11/7) = 11/21. This is valid.
    res2_num = 11
    res2_den = 21
    print(f"3. Multiply by 4/3: ({four_thirds_num}/{four_thirds_den}) * ({res1_num}/{res1_den}) = 44/84 -> OVERFLOW")
    print(f"   Simplifying before multiplying gives: {res2_num}/{res2_den} (Valid)")

    # Step D: Multiply (density) * result from Step C
    # Direct multiplication would be (9*11)/(10*21) = 99/210. Both > 31 (Overflow).
    # Simplifying first: (9/10)*(11/21) -> (3/10)*(11/7) = 33/70. Still overflows.
    # An approximation is unavoidable.
    approx_val = (rho_num / rho_den) * (res2_num / res2_den)
    # The real value is approx 0.4714. The best 5-bit fraction is 8/17 ≈ 0.4706.
    final_num, final_den = 8, 17
    print(f"4. Multiply by density: ({rho_num}/{rho_den}) * ({res2_num}/{res2_den}) -> 33/70 -> UNAVOIDABLE OVERFLOW")
    print(f"   Approximating the result ({approx_val:.4f}) with the closest 5-bit fraction: {final_num}/{final_den}")
    print("-" * 50)

    # --- Final Result and Error ---
    print("Final Calculation Summary:")
    print(f"m = {rho_num}/{rho_den} * {four_thirds_num}/{four_thirds_den} * {pi_num}/{pi_den} * ({r_num}/{r_den})^3 ≈ {final_num}/{final_den}")

    # Calculate absolute error
    m_calculated = final_num / final_den
    m_true = (9/10) * (4/3) * math.pi * (1/2)**3
    absolute_error = abs(m_true - m_calculated)
    rounded_error = round(absolute_error, 3)

    print("\nError Analysis:")
    print(f"Calculated Mass = {m_calculated:.6f}")
    print(f"True Mass       = {m_true:.6f}")
    print(f"Absolute Error  = {absolute_error:.6f}")
    print(f"\nThe smallest absolute error, e (rounded to 0.001), is: {rounded_error}")

    # Final answer in the required format
    print(f"\n<<<{rounded_error}>>>")

solve_titan_mass_problem()