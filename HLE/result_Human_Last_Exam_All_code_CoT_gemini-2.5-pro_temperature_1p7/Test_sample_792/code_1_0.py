import math

def main():
    """
    Solves the physics problem using the simulated Titan computer architecture.
    """

    print("--- Titan Computation Analysis ---")
    
    # 1. State the formula
    print("\nThe governing formula for the required force is F = (d * m * g) / h, where m is the rock's mass.")
    print("This must be calculated using Titan's 5-bit fractional arithmetic (all numerators/denominators <= 31).\n")

    # 2. Define constants and approximations
    print("Step 1: Define constants and select approximations.")
    d_n, d_d = 20, 1   # distance (aiming for center of 19-21m target)
    h_n, h_d = 10, 1   # height
    r_n, r_d = 1, 2    # radius (0.5 cm)
    rho_n, rho_d = 9, 10 # density (0.9 kg/cm^3)

    # Explain the choice of approximations for pi and g
    print("To keep all calculations within the 5-bit limit, we must use approximations that simplify well.")
    print("Using pi=3/1 and g=10/1 is necessary to prevent intermediate values from exceeding 31.")
    pi_n, pi_d = 3, 1
    g_n, g_d = 10, 1
    
    print(f"\nApproximations Used:")
    print(f"  d = {d_n}/{d_d} m")
    print(f"  h = {h_n}/{h_d} m")
    print(f"  pi = {pi_n}/{pi_d}")
    print(f"  g = {g_n}/{g_d} m/s^2\n")

    # 3. Calculate mass m
    print("Step 2: Calculate the mass (m) of the rock following Titan rules.")
    print("m = (4/3) * pi * r^3 * rho = (4/3) * (3/1) * (1/2)^3 * (9/10)")
    
    # (4/3) * (3/1) = 12/3 -> 4/1
    print("  a) (4/3) * (3/1) = 12/3. This simplifies to 4/1.")
    # r^3 = (1/2)^3 = 1/8
    print("  b) The radius cubed is (1/2)^3 = 1/8.")
    # m = (4/1) * (1/8) * (9/10)
    print("  c) The equation becomes m = (4/1) * (1/8) * (9/10).")
    # (4/1)*(1/8) = 4/8 -> 1/2
    print("  d) (4/1) * (1/8) = 4/8. This simplifies to 1/2.")
    # m = (1/2) * (9/10) = 9/20
    print("  e) (1/2) * (9/10) = 9/20.")
    m_n, m_d = 9, 20
    print(f"Resulting mass: m = {m_n}/{m_d} kg.\n")

    # 4. Calculate Force F
    print("Step 3: Calculate the required force (F).")
    print("F = (d/h) * m * g")
    
    # d/h = (20/1)/(10/1) = 2/1
    print(f"  a) (d/h) = ({d_n}/{d_d}) / ({h_n}/{h_d}) = 2/1.")
    d_over_h_n, d_over_h_d = 2, 1
    
    # F = (2/1) * (9/20) * (10/1)
    print(f"  b) The equation becomes F = ({d_over_h_n}/{d_over_h_d}) * ({m_n}/{m_d}) * ({g_n}/{g_d}).")
    
    # temp = (2/1) * (9/20) = 18/20 -> 9/10
    print(f"  c) ({d_over_h_n}/{d_over_h_d}) * ({m_n}/{m_d}) = 18/20. This simplifies to 9/10.")
    temp_n, temp_d = 9, 10
    
    # F = (9/10) * (10/1)
    print(f"  d) The final step is F = ({temp_n}/{temp_d}) * ({g_n}/{g_d}).")
    print("     A direct multiplication (9*10=90) would exceed the 31-limit.")
    print("     By rule 5, we apply algebraic simplification first, cancelling the common factor of 10.")
    print("     The calculation becomes (9/1) * (1/1) = 9/1.")
    F_n, F_d = 9, 1
    print(f"Resulting force: F = {F_n}/{F_d} N.\n")

    # 5. Show final equation
    print("--- Final Titan Equation ---")
    print(f"({d_n}/{d_d} m / {h_n}/{h_d} m) * ({m_n}/{m_d} kg) * ({g_n}/{g_d} m/s^2) = {F_n}/{F_d} N\n")

    # 6. Calculate Error
    print("--- Error Analysis ---")
    print("To find the error, we compare the Titan result to the high-precision 'true' value.")
    # High-precision values
    true_pi = math.pi
    true_g = 9.8
    r_val = 0.5
    rho_val = 0.9
    d_val = 20.0
    h_val = 10.0
    
    # True mass
    true_m = (4.0/3.0) * true_pi * (r_val**3) * rho_val
    # True force for d=20m
    true_F = (d_val/h_val) * true_m * true_g
    
    titan_F_val = F_n / F_d

    # Check if the calculated force hits the lion (target 19m-21m)
    F_min_true = (19.0/h_val) * true_m * true_g
    F_max_true = (21.0/h_val) * true_m * true_g
    
    print(f"True required force range to hit the lion (19m-21m): [{F_min_true:.5f}, {F_max_true:.5f}] N.")
    print(f"The Titan-calculated force is {titan_F_val:.1f} N, which is within this range.\nThus, the answer is Y[e].\n")

    print(f"True Force (for d=20m): {true_F:.5f} N")
    print(f"Titan Force:             {titan_F_val:.5f} N")
    
    error = abs(titan_F_val - true_F)
    print(f"The smallest absolute error (e) is |{titan_F_val:.5f} - {true_F:.5f}| = {error:.5f}")

if __name__ == "__main__":
    main()
