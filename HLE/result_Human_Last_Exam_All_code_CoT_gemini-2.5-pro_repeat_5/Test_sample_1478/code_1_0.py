import math

def find_best_fractional_approximation(value):
    """
    Finds the best fractional approximation for a given value that fits
    the Titan 6-bit architecture (numerator and denominator <= 63).
    The format is n/d * 10^e.
    """
    if value == 0:
        return 0, 1, 0, 0

    # Find the best exponent to bring the mantissa into a reasonable range
    # A good range for finding simple fractions is [1, 10)
    exponent = int(math.floor(math.log10(abs(value))))
    mantissa = value / (10**exponent)

    best_approx = {
        'n': 1,
        'd': 1,
        'exp': exponent,
        'error': float('inf')
    }

    # Search for the best fraction n/d for the mantissa
    for d in range(1, 64):
        n_float = mantissa * d
        # Test two candidate numerators: floor and ceiling
        for n_candidate in [math.floor(n_float), math.ceil(n_float)]:
            n = int(n_candidate)
            if n == 0 or n > 63:
                continue

            current_val = (n / d) * (10**exponent)
            error = abs(current_val - value) / abs(value)

            if error < best_approx['error']:
                best_approx['n'] = n
                best_approx['d'] = d
                best_approx['exp'] = exponent
                best_approx['error'] = error
                
    # Check if shifting the exponent by 1 yields a better approximation
    # for the new mantissa. e.g. 9.8 can be 98/10 * 10^0 or 98/1 * 10^-1
    mantissa_shifted = mantissa * 10
    exponent_shifted = exponent - 1
    for d in range(1, 64):
        n_float = mantissa_shifted * d
        for n_candidate in [math.floor(n_float), math.ceil(n_float)]:
            n = int(n_candidate)
            if n == 0 or n > 63:
                continue
            
            current_val = (n / d) * (10**exponent_shifted)
            error = abs(current_val - value) / abs(value)
            
            if error < best_approx['error']:
                best_approx['n'] = n
                best_approx['d'] = d
                best_approx['exp'] = exponent_shifted
                best_approx['error'] = error

    return best_approx['n'], best_approx['d'], best_approx['exp'], best_approx['error']

def solve_pandora_gravity():
    """
    Solves the physics problem based on the Titan architecture.
    """
    # 1. Define constants using high-precision floats for "true" value calculation
    G_true = 6.67430e-11
    PI_true = 3.1415926535
    
    # 2. Define planet and probe parameters
    R_planet = 2000 * 1000  # 2000 km in meters
    rho_planet = 1.2 * 1000  # 1.2 metric tons/m^3 in kg/m^3
    m_probe = 50.0  # kg
    d_from_horizon = 1000.0 # 1 km in meters

    # 3. Calculate "true" mass and force
    m_planet_true = rho_planet * (4.0/3.0) * PI_true * (R_planet**3)
    
    # Check if Schwarzschild radius is negligible
    # c_sq = (299792458.0)**2
    # rs_true = (2 * G_true * m_planet_true) / c_sq 
    # Since rs_true is very small (~0.04m), we can neglect it.
    r_total_true = d_from_horizon 
    
    F_true = G_true * m_probe * m_planet_true / (r_total_true**2)

    # 4. Simulate Titan calculation with fractional constants
    # Using the most accurate 6-bit friendly fractions
    # G = 20/3 * 10^-11
    G_n, G_d, G_exp = 20, 3, -11 
    # pi = 22/7
    PI_n, PI_d = 22, 7
    # rho = 1200 = 12/10 * 10^3
    rho_n, rho_d, rho_exp = 12, 10, 3
    # R = 2e6 = 2/1 * 10^6 => R^3 = 8/1 * 10^18
    R3_n, R3_d, R3_exp = 8, 1, 18
    # m1 = 50 = 50/1
    m1_n, m1_d = 50, 1
    # d = 1e3 => d^2 = 1e6
    d2_exp = 6
    
    # F = G * m1 * (rho * 4/3 * pi * R^3) / d^2
    # This is the calculation the 'RED' instruction would perform on the expression.
    F_titan_num = (G_n * m1_n * rho_n * 4 * PI_n * R3_n)
    F_titan_den = (G_d * m1_d * rho_d * 3 * PI_d * R3_d)
    F_titan_exp = G_exp + rho_exp + R3_exp - d2_exp
    
    F_titan_val = (F_titan_num / F_titan_den) * (10**F_titan_exp)

    # 5. Find the best Titan-representable fraction for the calculated force
    best_n, best_d, best_exp, approx_error = find_best_fractional_approximation(F_titan_val)
    
    # Final representation on Titan
    F_final_titan = (best_n / best_d) * (10**best_exp)

    # 6. Calculate the final relative error compared to the true value
    final_error = abs(F_final_titan - F_true) / F_true
    
    # 7. Print the results
    print("Titan Assembly Code:")
    print("--------------------")
    print("# This program calculates the gravitational force on a probe near Pandora.")
    print("# We determine the Schwarzschild radius is negligible and simplify r = 1km.")
    print("")
    print("# Load constants and parameters into registers")
    print("MOV AX, 20/3e-11   # Gravitational Constant (G)")
    print("MOV BX, 50/1       # Probe mass (m1)")
    print("MOV CX, 12/10e3    # Planet density (rho)")
    print("MOV DX, 22/7       # Pi")
    print("")
    print("# Calculate force F = G * m1 * (rho * 4/3 * pi * R^3) / d^2")
    print("MUL AX, BX         # AX = G * m1")
    print("MUL AX, CX         # AX = G * m1 * rho")
    print("MUL AX, 4/3        # AX = G * m1 * rho * 4/3")
    print("MUL AX, DX         # AX = G * m1 * rho * 4/3 * pi")
    print("MUL AX, 8/1e18     # AX = G * m1 * rho * 4/3 * pi * R^3")
    print("DIV AX, 1/1e6      # AX = F")
    print("")
    print("# Reduce the expression in AX to a single fractional value")
    print("RED AX")
    print("--------------------")
    print(f"\nThe calculation yields a force that can be best approximated on Titan as:")
    print(f"Final equation: F = {best_n}/{best_d} x 10^{best_exp} N")
    
    # The final answer format required by the prompt
    print(f"\n<<<Y[{final_error:.1%}]>>>")

solve_pandora_gravity()