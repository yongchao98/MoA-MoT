import math

# Helper function to find the greatest common divisor for fraction simplification
def gcd(a, b):
    return math.gcd(a, b)

# Global cache for simplification to speed up repeated lookups
simplification_cache = {}

def simplify_value(value, exp):
    """
    Finds the best 5-bit fractional representation for a given decimal value.
    It searches for a fraction n/d (n,d <= 31) and an exponent adjustment p
    such that (n/d) * 10^p is the best approximation of the original value.
    """
    if (value, exp) in simplification_cache:
        return simplification_cache[(value, exp)]

    best_approx = (None, None, None)
    min_error = float('inf')

    # Search for the best fraction by iterating through possible denominators
    for d_s in range(1, 32):
        # We want n_s / d_s approx value
        n_s = round(value * d_s)
        if 0 < n_s < 32:
            current_error = abs(value - n_s / d_s)
            if current_error < min_error:
                min_error = current_error
                # Store as (num, den, exponent)
                best_approx = (n_s, d_s, exp)

    simplification_cache[(value, exp)] = best_approx
    return best_approx


def simplify(n, d, e):
    """
    Simplifies a Titan number (n, d, e) if n or d are out of the 5-bit range.
    """
    if n <= 31 and d <= 31:
        common = gcd(n, d)
        return (n // common, d // common, e)
    
    val = n / d
    
    # We need to find n_s/d_s * 10^p that approximates val
    # This is equivalent to finding n_s/d_s that approximates val / 10^p
    # Let's try a few exponent shifts to find the best mantissa
    best_option = None
    min_error = float('inf')

    for p_shift in range(-2, 3):
        target_val = val / (10**p_shift)
        # Find best fractional approximation for target_val
        n_s, d_s, _ = simplify_value(target_val, 0)
        
        if n_s is not None:
            # Calculate the error of this representation
            approx_val = (n_s / d_s) * (10**p_shift)
            error = abs(val - approx_val)
            if error < min_error:
                min_error = error
                best_option = (n_s, d_s, e + p_shift)

    return best_option

# --- Titan Arithmetic Operations ---
def multiply(t1, t2):
    n1, d1, e1 = t1
    n2, d2, e2 = t2
    return simplify(n1 * n2, d1 * d2, e1 + e2)

def add(t1, t2):
    n1, d1, e1 = t1
    n2, d2, e2 = t2
    
    # Align exponents
    if e1 > e2:
        n2 *= 10**(e1 - e2)
        new_e = e1
    elif e2 > e1:
        n1 *= 10**(e2 - e1)
        new_e = e2
    else:
        new_e = e1
        
    # Add fractions
    new_n = n1 * d2 + n2 * d1
    new_d = d1 * d2
    return simplify(new_n, new_d, new_e)

def val(t):
    """Returns the decimal value of a Titan number."""
    return (t[0] / t[1]) * (10**t[2])

# --- Problem Calculation ---

# 1. Define constants in Titan format (n, d, e)
# Note: Approximations are chosen to be representable by 5-bit integers.
PI = (22, 7, 0)
FOUR_THIRDS = (4, 3, 0)
G = (20, 3, -11) # Approx 6.67e-11
M_PROBE = (5, 1, 1) # 50 kg
RHO_SHELL = (3, 1, 2) # 300 kg/m^3
RADIUS_A = (2, 1, 6) # 2e6 m

# 2. Calculate Gravitational Force (F_g)
# F_g ≈ (4/3)*pi * G * m_probe * rho_shell * a
# This simplified formula is derived from Fg = G*M*m/r^2 with M≈(4/3)pi*a^3*rho and r≈a
term1 = multiply(FOUR_THIRDS, PI) # 4/3 * pi -> 88/21 -> simplify to 21/5
term2 = multiply(term1, G) # * G -> (21/5)*(20/3)*10^-11 = 28*10^-11
term3 = multiply(term2, M_PROBE) # * 50 -> 140*10^-10 -> simplify to (7/5)*10^-8
term4 = multiply(term3, RHO_SHELL) # * 300 -> (21/5)*10^-6
f_g_titan = multiply(term4, RADIUS_A) # * 2e6 -> 42/5 -> simplify to 25/3

# Correct simplification path:
# 4/3 * pi = 88/21 ≈ 4.19 -> simplify to (21,5,0)
C1 = (21, 5, 0)
# C1 * G = (21/5)*(20/3)*10^-11 = 420/15 * 10^-11 = 28*10^-11 -> simplify to (28,1,-11)
T1 = simplify(21*20, 5*3, -11)
# T1 * M_PROBE = 28*10^-11 * 50 = 140*10^-10 -> simplify to (7,5,-8)
T2 = simplify(T1[0]*M_PROBE[0], T1[1]*M_PROBE[1], T1[2]+M_PROBE[2])
# T2 * RHO_SHELL = (7/5)*10^-8 * 300 = 21/5 * 10^-6 -> simplify to (21,5,-6)
T3 = simplify(T2[0]*RHO_SHELL[0], T2[1]*RHO_SHELL[1], T2[2]+RHO_SHELL[2])
# T3 * RADIUS_A = (21/5)*10^-6 * 2e6 = 42/5 -> simplify to (25,3,0)
f_g_titan = simplify(T3[0]*RADIUS_A[0], T3[1]*RADIUS_A[1], T3[2]+RADIUS_A[2])


# 3. Calculate Deceleration Force (F_decel)
# a = v_i^2 / (2*d) = 300^2 / (2*5000) = 9 m/s^2
# F_decel = m * a = 50 * 9 = 450 N
# Represent 450 N in Titan format: 4.5 * 10^2 = (9/2) * 10^2
f_decel_titan = (9, 2, 2)

# 4. Calculate Total Rocket Force (F_rocket)
# To add, we must align exponents.
# F_g = 25/3 = 8.333 = 0.08333 * 10^2. Best Titan fraction for 0.08333 is 1/12.
f_g_aligned = (1, 12, 2)
# F_rocket = (1/12 + 9/2) * 10^2 = ( (2+108)/24 ) * 10^2 = 110/24 * 10^2
# Simplify 110/24 * 10^2. 110/24 = 11/2.4 = 1.1/0.24.
# 110/24 = 4.5833... = 4.5833 * 10^0. This can be written as 0.45833 * 10^1.
# Best Titan fraction for 0.45833 is 11/24.
# So we have (11/24) * 10^1 * 10^2 = (11/24) * 10^3
f_rocket_titan = simplify(110, 24, 2)


# 5. Print the final equation
print("[Precise Landing with Superconducting Computer]")
print("\nFinal force calculation based on Titan's 5-bit fractional arithmetic:")
print("\nTarget Equation: F_rocket = F_deceleration + F_gravity")

# Print the components with their values and Titan representations
print(f"\n1. Deceleration Force (F_decel):")
print(f"   - Titan Representation: ({f_decel_titan[0]}/{f_decel_titan[1]}) * 10^{f_decel_titan[2]}")
print(f"   - Value: {val(f_decel_titan):.4f} N")

print(f"\n2. Gravitational Force (F_g):")
print(f"   - Titan Representation: ({f_g_titan[0]}/{f_g_titan[1]}) * 10^{f_g_titan[2]}")
print(f"   - Value: {val(f_g_titan):.4f} N")

print(f"\n3. Total Required Rocket Force (F_rocket):")
print(f"   - Titan Representation: ({f_rocket_titan[0]}/{f_rocket_titan[1]}) * 10^{f_rocket_titan[2]}")
print(f"   - Value: {val(f_rocket_titan):.4f} N")

# Print the final equation with calculated values
print("\n--- FINAL EQUATION ---")
print(f"{val(f_rocket_titan):.4f} N = {val(f_decel_titan):.4f} N + {val(f_g_titan):.4f} N")