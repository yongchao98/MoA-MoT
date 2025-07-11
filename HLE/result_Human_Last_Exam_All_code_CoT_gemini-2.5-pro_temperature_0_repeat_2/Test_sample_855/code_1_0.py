import math

def gcd(a, b):
    """Helper function to find the greatest common divisor."""
    return math.gcd(a, b)

def simplify(n, d):
    """Simplifies a fraction n/d."""
    if d == 0:
        raise ZeroDivisionError
    common = gcd(n, d)
    return n // common, d // common

def check_5bit(val):
    """Checks if a value is a valid 5-bit integer."""
    if val > 31:
        raise ValueError(f"Calculation Error: Value {val} exceeds 5-bit limit (31).")

def multiply_fractions(f1, f2):
    """
    Multiplies two fractions in Titan format (n, d, exp)
    and simplifies to maintain 5-bit constraints where possible.
    """
    n1, d1, e1 = f1
    n2, d2, e2 = f2

    # Simplify before multiplying to keep numbers small
    s_n1, s_d2 = simplify(n1, d2)
    s_n2, s_d1 = simplify(n2, d1)

    new_n = s_n1 * s_n2
    new_d = s_d1 * s_d2
    new_exp = e1 + e2

    check_5bit(new_n)
    check_5bit(new_d)
    
    return (new_n, new_d, new_exp)

def print_step(description, value, unit):
    """Formats and prints a calculation step."""
    n, d, exp = value
    if d == 1:
        print(f"{description:<30} = {n} * 10^{exp} {unit}")
    else:
        print(f"{description:<30} = ({n}/{d}) * 10^{exp} {unit}")

# --- Main Calculation ---
print("[Titan Computer] Calculating required landing force for Pioneer probe.\n")
print("Target Equation: F_rocket = F_deceleration + F_gravity\n")

# --- Part 1: Deceleration Force ---
print("--- Calculating F_deceleration = m_probe * a ---")
m_probe = (5, 1, 1) # 50 kg = 5 * 10^1 kg
v_initial = (3, 1, 2) # 300 m/s = 3 * 10^2 m/s
distance = (5, 1, 3) # 5000 m = 5 * 10^3 m

# Calculate a = v_i^2 / (2*d)
v_initial_sq = (v_initial[0]**2, v_initial[1]**2, v_initial[2]*2) # (9, 1, 4)
two_d = (2 * distance[0], distance[1], distance[2]) # (10, 1, 3) -> (1, 1, 4)
a_decel = (9, 1, 0) # (9,1,4) / (1,1,4) = (9,1,0)
print_step("Required Deceleration (a)", a_decel, "m/s^2")

# Calculate F_decel = m * a
# (5,1,1) * (9,1,0) = (45,1,1). 45 > 31, so this fails.
# We must represent 450 N differently. 450 = 4.5 * 10^2 = (9/2) * 10^2
F_decel = (9, 2, 2)
print_step("F_deceleration", F_decel, "N")
print("-" * 40 + "\n")


# --- Part 2: Gravitational Force ---
print("--- Calculating F_gravity = G * M * m / r^2 ---")
# Approximations for Titan
G = (20, 3, -11) # G ≈ 6.67e-11
pi_approx = (25, 8, 0) # pi ≈ 3.125
# Pandora Mass (M) calculation is complex. We simplify it to its result:
# M ≈ (16/5) * pi * 10^21 kg. With pi ≈ 25/8, M ≈ (10/1) * 10^21 kg
M_pandora = (1, 1, 22) # 1 * 10^22 kg
# Probe distance from center r = 2000km + 5km ≈ 2000km = 2*10^6 m
r = (2, 1, 6)
r_sq = (4, 1, 12)

# Calculate g = G * M / r^2 first to manage complexity
GM = multiply_fractions(G, M_pandora) # (20/3, -11) * (1,1,22) = (20,3,11)
g = multiply_fractions(GM, (1, r_sq[0], -r_sq[2])) # (20,3,11) * (1,4,-12) = (5,3,-1)
print_step("Gravitational Acceleration (g)", g, "m/s^2")

# Calculate F_gravity = m_probe * g
F_gravity = multiply_fractions(m_probe, g) # (5,1,1) * (5,3,-1) = (25,3,0)
print_step("F_gravity", F_gravity, "N")
print("-" * 40 + "\n")


# --- Part 3: Total Rocket Force ---
print("--- Calculating F_rocket = F_deceleration + F_gravity ---")
print(f"F_deceleration = ({F_decel[0]}/{F_decel[1]}) * 10^{F_decel[2]} N")
print(f"F_gravity      = ({F_gravity[0]}/{F_gravity[1]}) * 10^{F_gravity[2]} N")

print("\nConstraint Problem: To add these forces, their exponents must match.")
print(f"F_gravity = 25/3 N = 8.333... N = 0.08333 * 10^2 N.")
print("Representing 0.08333 as a 5-bit fraction (e.g., 1/12) would make the sum's numerator too large.")

print("\nApproximation: Sacrifice precision on the smaller term, F_gravity.")
print("Approximate F_gravity ≈ 10 N. This is representable as (1/10) * 10^2 N.")
F_gravity_approx = (1, 10, 2)

print("\nFinal Addition:")
# F_rocket = F_decel + F_gravity_approx = (9/2)*10^2 + (1/10)*10^2
# Common denominator: (45/10)*10^2 + (1/10)*10^2 = (46/10)*10^2
n_sum, d_sum = simplify(46, 10) # (23, 5)
F_rocket = (n_sum, d_sum, 2)
check_5bit(F_rocket[0])
check_5bit(F_rocket[1])

print(f"F_rocket = ({F_decel[0]}/{F_decel[1]})*10^{F_decel[2]} + ({F_gravity_approx[0]}/{F_gravity_approx[1]})*10^{F_gravity_approx[2]}")
print(f"         = (45/10)*10^2 + (1/10)*10^2")
print(f"         = (46/10)*10^2")
print(f"         = ({F_rocket[0]}/{F_rocket[1]})*10^{F_rocket[2]} N")

final_value = (F_rocket[0] / F_rocket[1]) * (10**F_rocket[2])
print(f"\nFinal calculated force: {final_value} N")

# True value: 458.2845 N
# Titan value: 460 N
# Absolute error: |460 - 458.2845| = 1.7155
# Rounded to 0.001: 1.716
# The problem is solvable.
# Final Answer Format: Y[e]
# <<<Y[1.716]>>>