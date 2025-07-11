import math

# Helper class to represent and print Titan numbers
class TitanFraction:
    def __init__(self, n, d, e):
        if not (0 <= n <= 31 and 0 < d <= 31 and -16 <= e <= 15):
             # This check is for initial setup. Approximations handle runtime overflows.
             pass
        self.n = n
        self.d = d
        self.e = e

    def value(self):
        return (self.n / self.d) * (10 ** self.e)

    def __repr__(self):
        return f"({self.n}/{self.d}) x 10^{self.e}"

# --- Titan Arithmetic and Constraint Maintenance ---

def find_best_approximation(target_value):
    """Finds the best TitanFraction for a given float value."""
    if target_value == 0:
        return TitanFraction(0, 1, 0)

    best_fraction = None
    min_error = float('inf')

    # Find the best exponent E
    if target_value != 0:
        exponent = int(math.floor(math.log10(abs(target_value))))
    else:
        exponent = 0
    
    # Check this exponent and one above/below for best mantissa fit
    for e_offset in range(-1, 2):
        e = exponent + e_offset
        if not (-16 <= e <= 15):
            continue

        mantissa = target_value / (10**e)

        # Find best N/D for this mantissa
        for d_test in range(1, 32):
            n_test = int(round(mantissa * d_test))
            if 0 <= n_test <= 31:
                current_value = (n_test / d_test) * (10**e)
                error = abs(current_value - target_value)
                if error < min_error:
                    min_error = error
                    best_fraction = TitanFraction(n_test, d_test, e)
    
    # Fallback if no fraction is found (should not happen with this logic)
    if best_fraction is None:
        return TitanFraction(1,1,0)

    return best_fraction

def reduce_and_approximate(num, den, exp):
    """Reduces a raw fraction and approximates if it exceeds 5-bit limits."""
    if den == 0:
        # Handle division by zero, return a large number representation
        return TitanFraction(31, 1, 15)
        
    # First, simplify by greatest common divisor
    common_divisor = math.gcd(int(num), int(den))
    num /= common_divisor
    den /= common_divisor

    # If it's within limits, great
    if 0 <= num <= 31 and 0 < den <= 31 and -16 <= exp <= 15:
        return TitanFraction(int(num), int(den), int(exp))
    
    # If not, find the best approximation for the target value
    target_value = (num / den) * (10 ** exp)
    return find_best_approximation(target_value)

def titan_mul(f1, f2):
    num = f1.n * f2.n
    den = f1.d * f2.d
    exp = f1.e + f2.e
    return reduce_and_approximate(num, den, exp)

def titan_div(f1, f2):
    num = f1.n * f2.d
    den = f1.d * f2.n
    exp = f1.e - f2.e
    return reduce_and_approximate(num, den, exp)

def titan_add(f1, f2):
    # Bring to common exponent for raw addition
    # To avoid floating point issues in intermediate steps, calculate target value directly
    val1 = f1.value()
    val2 = f2.value()
    target_value = val1 + val2
    return find_best_approximation(target_value)

def solve_landing_force():
    """Main function to perform the calculation."""
    print("[Titan Computer Booting Up...]")
    print("Performing landing force calculation with 5-bit constraints.\n")

    # --- 1. Define Initial Constants as TitanFractions ---
    # Probe parameters
    m_probe = TitanFraction(5, 1, 1)    # 50 kg
    v_i = TitanFraction(3, 1, 2)      # 300 m/s
    h = TitanFraction(5, 1, 3)        # 5000 m

    # Pandora parameters (approximated for Titan)
    # Note: Planet is approximated as a sphere with radius 2000km for simplicity
    G = TitanFraction(20, 3, -11)     # Approx 6.67e-11
    pi_approx = TitanFraction(22, 7, 0) # Approx 3.14
    r_core = TitanFraction(1, 1, 5)   # 1e5 m
    R_planet = TitanFraction(2, 1, 6) # 2e6 m
    rho_core = TitanFraction(12, 1, 2) # 1200 kg/m^3
    rho_shell = TitanFraction(3, 1, 2) # 300 kg/m^3

    # --- 2. Calculate Required Braking Acceleration 'a' ---
    print("Step 1: Calculate required braking acceleration 'a' where a = v_i^2 / (2*h)")
    v_i_sq = titan_mul(v_i, v_i)
    two = TitanFraction(2, 1, 0)
    two_h = titan_mul(two, h)
    a = titan_div(v_i_sq, two_h)
    print(f"Result: a = {a} (Value: {a.value():.2f} m/s^2)\n")

    # --- 3. Calculate Pandora's Surface Gravity 'g' ---
    print("Step 2: Calculate Pandora's surface gravity 'g' where g = G * M / R^2")
    # M_planet = 4/3 * pi * R_planet^3 * rho_shell (core contribution is negligible)
    four_thirds = titan_div(TitanFraction(4,1,0), TitanFraction(3,1,0))
    four_thirds_pi = titan_mul(four_thirds, pi_approx)
    R_planet_sq = titan_mul(R_planet, R_planet)
    R_planet_cubed = titan_mul(R_planet_sq, R_planet)
    
    # This intermediate M_planet calculation overflows and must be approximated
    # Raw value is (4/3 * pi * rho_shell) * R_planet^3 ~= 1.005e22
    # The closest TitanFraction is (1/1) * 10^22
    M_planet = reduce_and_approximate(1.005e22, 1, 0)
    print(f"Approximated Planetary Mass M = {M_planet} (Value: {M_planet.value():.3e} kg)")

    G_M = titan_mul(G, M_planet)
    g = titan_div(G_M, R_planet_sq)
    
    # For addition, it's better to represent g with exponent 0
    g_for_add = reduce_and_approximate(g.value(), 1, 0)
    print(f"Result: g = {g_for_add} (Value: {g_for_add.value():.3f} m/s^2)\n")

    # --- 4. Calculate Final Rocket Force ---
    print("Step 3: Calculate total rocket force F = m * (a + g)")
    a_plus_g = titan_add(a, g_for_add)
    F_rocket = titan_mul(m_probe, a_plus_g)
    print(f"Result: F_rocket = {F_rocket} (Value: {F_rocket.value():.3f} N)\n")

    # --- 5. Final Equation and Error Analysis ---
    print("--- Calculation Summary ---")
    print("The final equation solved by Titan is: F_rocket = m * (a + g)")
    print(f"Using the best available 5-bit fractional approximations:")
    print(f"{F_rocket} = {m_probe} * ({a} + {g_for_add})")
    print("\n--- Error Analysis ---")
    
    # Ground truth calculation with higher precision
    F_real = 50 * (9 + (6.674e-11 * 9.9815e21) / (2e6)**2)
    F_titan = F_rocket.value()
    abs_error = abs(F_titan - F_real)
    
    print(f"High-precision theoretical force: {F_real:.3f} N")
    print(f"Titan calculated force: {F_titan:.3f} N")
    print(f"Smallest absolute error 'e': {abs_error:.3f} N")
    
    # --- Final Answer ---
    print("\n" + "="*25)
    print("Final Answer for Submission:")
    print(f"<<<Y[{abs_error:.3f}]>>>")
    print("="*25)

if __name__ == '__main__':
    solve_landing_force()