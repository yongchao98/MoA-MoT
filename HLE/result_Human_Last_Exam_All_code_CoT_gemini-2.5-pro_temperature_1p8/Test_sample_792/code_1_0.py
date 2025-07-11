import math

# This simulation solves the physics problem using the computational rules of
# the Titan supercomputer. All calculations are performed with fractional
# arithmetic where numerators and denominators are 5-bit integers (0-31).

# --- Helper function for printing fractions ---
def print_fraction(n, d):
    """Prints a fraction in n/d format."""
    return f"{int(n)}/{int(d)}"

# --- Main Calculation ---
def titan_force_calculation():
    """
    Calculates the required force using Titan's 5-bit fractional arithmetic.
    """
    # Step 1: Define constants and choose safe approximations
    # Physical quantities are represented as fractions.
    # r = 0.5 cm -> 1/2 cm
    r_n, r_d = 1, 2
    # rho = 0.9 kg/cm^3 -> 9/10 kg/cm^3
    rho_n, rho_d = 9, 10

    # Approximations for physical constants are chosen to keep calculations
    # within 5-bit integer limits (0-31).
    # pi is approximated as 3/1.
    pi_n, pi_d = 3, 1
    # g is approximated as 10/1 (from true value of ~9.8).
    g_n, g_d = 10, 1

    print("--- Titan Calculation for Required Force F ---")
    print("The formula for the force is: F = 2 * m * g")
    print(f"We will use the approximations: pi = {print_fraction(pi_n, pi_d)}, g = {print_fraction(g_n, g_d)}\n")

    # Step 2: Calculate the mass 'm' of the rock
    # m = Volume * density = (4/3 * pi * r^3) * rho
    print("Step A: Calculating Mass (m)")

    # 1. Calculate Volume component V_part1 = 4/3 * pi
    vp1_n_res = 4 * pi_n
    vp1_d_res = 3 * pi_d
    vp1_n_s, vp1_d_s = 4, 1 # Simplified from 12/3
    print(f"  1. (4/3) * {print_fraction(pi_n, pi_d)} = {print_fraction(vp1_n_res, vp1_d_res)} => simplify to {print_fraction(vp1_n_s, vp1_d_s)}")

    # 2. Calculate r^3
    r3_n, propensity = r_n**3, r_d**3
    print(f"  2. r^3 = ({print_fraction(r_n, r_d)})^3 = {print_fraction(r3_n, r3_d)}")

    # 3. Calculate full Volume V = V_part1 * r^3
    v_n_res = vp1_n_s * r3_n
    v_d_res = vp1_d_s * r3_d
    v_n_s, v_d_s = 1, 2 # Simplified from 4/8
    print(f"  3. Volume V = {print_fraction(vp1_n_s, vp1_d_s)} * {print_fraction(r3_n, r3_d)} = {print_fraction(v_n_res, v_d_res)} => simplify to {print_fraction(v_n_s, v_d_s)}")

    # 4. Calculate Mass m = V * rho
    m_n, m_d = v_n_s * rho_n, v_d_s * rho_d
    print(f"  4. Mass m = {print_fraction(v_n_s, v_d_s)} * {print_fraction(rho_n, rho_d)} = {print_fraction(m_n, m_d)}")
    print(f"Result: Calculated Mass m = {print_fraction(m_n, m_d)} kg.\n")

    # Step 3: Calculate the Force 'F'
    # F = 2 * m * g
    print("Step B: Calculating Force (F)")

    # 1. Calculate 2 * m
    fm_n_res = 2 * m_n
    fm_d_res = 1 * m_d
    fm_n_s, fm_d_s = 9, 10 # Simplified from 18/20
    print(f"  1. (2/1) * {print_fraction(m_n, m_d)} = {print_fraction(fm_n_res, fm_d_res)} => simplify to {print_fraction(fm_n_s, fm_d_s)}")

    # 2. Calculate final Force F = (2*m) * g
    # We simplify by cancelling common factors to prevent overflow. (9/10) * (10/1) = 9/1
    F_n, F_d = 9, 1
    print(f"  2. ({print_fraction(fm_n_s, fm_d_s)}) * {print_fraction(g_n, g_d)} => simplify to {print_fraction(F_n, F_d)}")
    
    print("\n--- Final Equation and Result ---")
    print(f"F = (2/1) * ({m_n}/{m_d}) * ({g_n}/{g_d}) = {print_fraction(F_n, F_d)} N")

# Run the calculation
titan_force_calculation()

# Final Answer calculation based on problem format Y[e]
# True force = 2 * (0.15 * pi) * 9.8
true_force = 2 * (0.15 * math.pi) * 9.8
titan_force = 9
absolute_error = abs(true_force - titan_force)
# print(f"\n<<<Y[{absolute_error:.3f}]>>>") # This is for generating the final response line