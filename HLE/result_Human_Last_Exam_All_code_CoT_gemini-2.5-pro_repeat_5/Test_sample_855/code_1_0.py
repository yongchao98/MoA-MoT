import math

# [Precise Landing with Superconducting Computer]
# This script simulates the Titan computer to calculate the required landing force.

# ==============================================================================
# Helper functions for Titan 5-bit Fractional Arithmetic
# ==============================================================================

def gcd(a, b):
    """Computes the greatest common divisor of two integers."""
    while b:
        a, b = b, a % b
    return abs(a)

def find_best_fraction(value, max_den=31):
    """
    Finds the best fractional approximation p/q for a given decimal value,
    constrained by Titan's 5-bit register size (p, q <= 31).
    """
    if value == 0:
        return 0, 1
    best_frac = (1, 1)
    min_error = float('inf')
    for den in range(1, max_den + 1):
        num = int(round(value * den))
        if 0 < num <= 31:
            error = abs(value - num / den)
            if error < min_error:
                min_error = error
                best_frac = (num, den)
    return best_frac

def normalize(num, den, exp):
    """
    Normalizes a number (num/den * 10^exp) to fit Titan's constraints.
    It adjusts the fraction and exponent while keeping num and den within the
    5-bit range [0, 31] and the exponent in the signed 5-bit range [-16, 15].
    """
    if num == 0:
        return 0, 1, 0

    common = gcd(num, den)
    num //= common
    den //= common

    if num <= 31 and den <= 31:
        return num, den, exp

    value = num / den
    k = 0
    if value > 0:
        k = math.floor(math.log10(value))

    target_val_for_frac = value / (10**k)
    new_exp = exp + k
    p, q = find_best_fraction(target_val_for_frac, 31)

    if not -16 <= new_exp <= 15:
        print(f"Warning: Exponent {new_exp} is out of [-16, 15] range.")

    return p, q, new_exp

def titan_mul(t1, t2):
    """Multiplies two Titan numbers t1 * t2."""
    n1, d1, e1 = t1
    n2, d2, e2 = t2
    new_num = n1 * n2
    new_den = d1 * d2
    new_exp = e1 + e2
    return normalize(new_num, new_den, new_exp)

def titan_add(t1, t2):
    """Adds two Titan numbers t1 + t2."""
    n1, d1, e1 = t1
    n2, d2, e2 = t2

    if e1 > e2:
        diff = e1 - e2
        d2_new = d2 * (10**diff)
        n2_new = n2
        final_exp = e1
        num = n1 * d2_new + n2_new * d1
        den = d1 * d2_new
    elif e2 > e1:
        diff = e2 - e1
        d1_new = d1 * (10**diff)
        n1_new = n1
        final_exp = e2
        num = n1_new * d2 + n2 * d1_new
        den = d2 * d1_new
    else:
        final_exp = e1
        num = n1 * d2 + n2 * d1
        den = d1 * d2

    return normalize(num, den, final_exp)

# ==============================================================================
# Main Calculation
# ==============================================================================

# Define constants in Titan's fractional format (num, den, exp)
G = (20, 3, -11)          # Approx 6.67e-11 m^3 kg^-1 s^-2
FOUR_THIRDS = (4, 3, 0)
PI = (22, 7, 0)           # Approx 3.142
R_PLANET = (2, 1, 6)      # Pandora's radius ~2,000,000 m
RHO_SHELL = (3, 1, 2)     # Shell density 300 kg/m^3
M_PROBE = (5, 1, 1)       # Probe mass 50 kg
A_DECEL = (9, 1, 0)       # Required deceleration 9 m/s^2

# --- Step 1: Calculate Pandora's surface gravity, g ---
# g = G * (4/3) * pi * R_planet * rho_shell
g_term1 = titan_mul(G, FOUR_THIRDS)
g_term2 = titan_mul(g_term1, PI)
g_term3 = titan_mul(g_term2, R_PLANET)
g_pandora = titan_mul(g_term3, RHO_SHELL)

# --- Step 2: Calculate the total required acceleration (g + a) ---
total_accel = titan_add(g_pandora, A_DECEL)

# --- Step 3: Calculate the final rocket force F = m * (g + a) ---
f_rocket = titan_mul(M_PROBE, total_accel)

# --- Step 4: Output the equation with final Titan numbers ---
print("--- Titan Calculation Steps ---")
print("Equation: F_rocket = m_probe * (g_pandora + a_decel)")
print("where g_pandora = G * (4/3) * pi * R_planet * rho_shell\n")

g_n, g_d, g_e = g_pandora
a_n, a_d, a_e = A_DECEL
f_n, f_d, f_e = f_rocket
m_n, m_d, m_e = M_PROBE

print("Final equation using Titan-normalized fractions:")
print(f"F_rocket = ({m_n}/{m_d})*10^{m_e} * ( ({g_n}/{g_d})*10^{g_e} + ({a_n}/{a_d})*10^{a_e} )")
print(f"F_rocket = ({f_n}/{f_d})*10^{f_e}")

# --- Step 5: Calculate the final numerical value and compare with a precise calculation ---
titan_force_value = (f_n / f_d) * (10**f_e)

# Precise calculation (using floating point numbers for comparison)
g_val = 6.674e-11 * (4/3) * math.pi * 2e6 * 300
real_force_value = 50 * (g_val + 9)
abs_error = abs(titan_force_value - real_force_value)

print(f"\n--- Results ---")
print(f"Final calculated force from Titan: {titan_force_value:.3f} N")
print(f"Precise force value (for comparison): {real_force_value:.3f} N")
print(f"Smallest absolute error: {abs_error:.3f}")

print(f"\n<<<Y[{abs_error:.3f}]>>>")