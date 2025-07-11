import math

def check_5bit(n, d):
    """Helper to check if a fraction is valid in Titan architecture."""
    if n > 31 or d > 31:
        raise ValueError(f"Constraint violation: {n}/{d}. Numerator or denominator exceeds 31.")

def simplify(n, d):
    """Helper to simplify a fraction using GCD."""
    if d == 0:
        raise ValueError("Denominator cannot be zero.")
    common_divisor = math.gcd(n, d)
    return n // common_divisor, d // common_divisor

def main():
    """
    Solves the Curious Monkey problem using simulated Titan 5-bit fractional arithmetic.
    """
    print("--- Titan Computation Plan ---")
    print("Goal: Calculate the force F to hit a target at x = 20m from a height of h = 10m.")
    print("Formula: F = (x / h) * m * g\n")

    # --- Step 1: Define Constants as Titan Fractions ---
    print("--- Step 1: Define Constants and Approximations ---")
    h_n, h_d = 10, 1
    x_n, x_d = 20, 1
    rho_n, rho_d = 9, 10 # 0.9 kg/cm^3
    r_n, r_d = 1, 2      # 0.5 cm
    four_thirds_n, four_thirds_d = 4, 3
    pi_approx_n, pi_approx_d = 3, 1     # Approximation for pi
    g_approx_n, g_approx_d = 10, 1    # Approximation for g (m/s^2)
    
    print(f"h = {h_n}/{h_d} m")
    print(f"x = {x_n}/{x_d} m")
    print(f"rho = {rho_n}/{rho_d} kg/cm^3")
    print(f"r = {r_n}/{r_d} cm")
    print(f"pi (approx) = {pi_approx_n}/{pi_approx_d}")
    print(f"g (approx) = {g_approx_n}/{g_approx_d} m/s^2\n")

    # --- Step 2: Calculate Mass (m) step-by-step ---
    print("--- Step 2: Calculate Mass m = rho * (4/3) * pi * r^3 ---")
    
    # r^3
    r3_n, r3_d = r_n**3, r_d**3
    print(f"r^3 = ({r_n}/{r_d})^3 = {r3_n}/{r3_d}")
    check_5bit(r3_n, r3_d)

    # Perform multiplication in an order that avoids overflow
    # Term A = rho * pi_approx
    print(f"Calculating (rho * pi): ({rho_n}/{rho_d}) * ({pi_approx_n}/{pi_approx_d})")
    a_n, a_d = rho_n * pi_approx_n, rho_d * pi_approx_d
    a_n, a_d = simplify(a_n, a_d)
    print(f"Result A = {a_n}/{a_d}")
    check_5bit(a_n, a_d)
    
    # Term B = (4/3) * r^3
    print(f"Calculating (4/3 * r^3): ({four_thirds_n}/{four_thirds_d}) * ({r3_n}/{r3_d})")
    b_n, b_d = four_thirds_n * r3_n, four_thirds_d * r3_d
    b_n, b_d = simplify(b_n, b_d)
    print(f"Result B = {b_n}/{b_d}")
    check_5bit(b_n, b_d)
    
    # m = Term A * Term B
    print(f"Calculating m = A * B: ({a_n}/{a_d}) * ({b_n}/{b_d})")
    m_n, m_d = a_n * b_n, a_d * b_d
    m_n, m_d = simplify(m_n, m_d)
    print(f"Final Mass m = {m_n}/{m_d} kg\n")
    check_5bit(m_n, m_d)

    # --- Step 3: Calculate Force (F) step-by-step ---
    print("--- Step 3: Calculate Force F = (x / h) * m * g ---")
    
    # Term C = x / h
    print(f"Calculating (x / h): ({x_n}/{x_d}) / ({h_n}/{h_d})")
    c_n, c_d = x_n * h_d, x_d * h_n
    c_n, c_d = simplify(c_n, c_d)
    print(f"Result C = {c_n}/{c_d}")
    check_5bit(c_n, c_d)

    # Term D = C * m
    print(f"Calculating C * m: ({c_n}/{c_d}) * ({m_n}/{m_d})")
    d_n, d_d = c_n * m_n, c_d * m_d
    d_n, d_d = simplify(d_n, d_d)
    print(f"Result D = {d_n}/{d_d}")
    check_5bit(d_n, d_d)

    # F = D * g
    print(f"Calculating F = D * g: ({d_n}/{d_d}) * ({g_approx_n}/{g_approx_d})")
    f_n, f_d = d_n * g_approx_n, d_d * g_approx_d
    f_n, f_d = simplify(f_n, f_d)
    f_approx = f_n / f_d
    print(f"Final Approximated Force F = {f_n}/{f_d} = {f_approx:.1f} N\n")
    check_5bit(f_n, f_d)
    
    # --- Step 4: Display Final Equation and Calculate Error ---
    print("--- Step 4: Final Equation and Error Analysis ---")
    print("The final equation with our fractional numbers is:")
    print(f"F = (x/h) * m * g")
    print(f"{f_n}/{f_d} = (({x_n}/{x_d}) / ({h_n}/{h_d})) * ({m_n}/{m_d}) * ({g_approx_n}/{g_approx_d})")
    print(f"{f_n}/{f_d} = ({c_n}/{c_d}) * ({m_n}/{m_d}) * ({g_approx_n}/{g_approx_d})\n")

    # Calculate true force for error comparison
    m_true = (0.9 * math.pi / 6)
    g_true = 9.8
    f_true = (20.0 / 10.0) * m_true * g_true
    error = abs(f_approx - f_true)
    
    print(f"True mass (high precision)      ~ {m_true:.5f} kg")
    print(f"True force (high precision)     ~ {f_true:.5f} N")
    print(f"Approximated force              = {f_approx:.5f} N")
    print(f"Absolute Error |F_true - F_approx| = {error:.5f} N\n")
    
    print("Conclusion: The calculation is possible on the Titan architecture.")
    print("The smallest absolute error we can produce with a valid calculation is {:.3f}.".format(error))
    
    # Final answer format
    print("\n<<<Y[{:.3f}]>>>".format(error))

if __name__ == '__main__':
    main()