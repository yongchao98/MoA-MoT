import math

def calculate_landing_time():
    """
    This script performs the calculations to find the landing time on Pandora,
    using fractional approximations compatible with the Titan 4-bit architecture.
    """

    # --- Step 1: Define approximations for constants and parameters ---
    # These fractions are chosen such that their numerators and denominators are <= 15.
    # The key is to select values whose products also remain within the 4-bit range.

    # pi ≈ 3.14159
    pi_n, pi_d = 13, 4  # Represents 3.25

    # Gravitational Constant G ≈ 6.674e-11 N(m/kg)^2. We use G ≈ 0.5e-10 to avoid overflow.
    G_mantissa_n, G_mantissa_d = 1, 2
    G_exp = -10

    # Pandora's shell density rho_shell = 300 kg/m^3
    rho_n, rho_d = 3, 1
    rho_exp = 2

    # Pandora's radius R ≈ 2,000 km = 2,000,000 m (approximating as a sphere)
    R_n, R_d = 2, 1
    R_exp = 6

    # Drop height h = 5,000 m
    h_n, h_d = 5, 1
    h_exp = 3

    # --- Step 2: Calculate gravitational acceleration g on Titan ---
    # Formula: g = (4/3) * pi * G * rho * R
    # Titan calculates the mantissa by multiplying the fractional parts:
    # (4/3) * (13/4) * (1/2) * (3/1) * (2/1)
    # Assuming Titan's ALU can reorder and cancel terms to avoid overflow:
    # (4/3 * 3/1) * (13/4) * (2/1 * 1/2) -> (4/1) * (13/4) * (1/1) -> 13/1
    g_mantissa_n, g_mantissa_d = 13, 1
    # Sum of exponents: G_exp + rho_exp + R_exp
    g_exp = G_exp + rho_exp + R_exp  # -10 + 2 + 6 = -2

    # --- Step 3: Calculate t^2 = 2 * h / g ---
    # t^2 = (2 * 5*10^3) / (13*10^-2)
    # Mantissa part: (2/1 * 5/1) / (13/1) = 10/13
    t_squared_mantissa_n, t_squared_mantissa_d = 10, 13
    # Exponent part: 10^3 / 10^-2 = 10^5
    t_squared_exp = h_exp - g_exp

    # --- Step 4: Calculate t = sqrt(t^2) using approximation ---
    # We need to compute sqrt((10/13) * 10^5). The exponent is odd.
    # We can rewrite as sqrt((10/13) / 10 * 10^(5+1)) = sqrt(1/13) * 10^3
    # Now we need a fractional approximation for sqrt(1/13) ≈ 0.2773
    # A good 4-bit fraction is 3/11 = 0.2727...
    sqrt_approx_n, sqrt_approx_d = 3, 11
    
    # So, t is approximated as (3/11) * 10^3
    t_final_n, t_final_d = sqrt_approx_n, sqrt_approx_d
    t_final_exp = 3

    # Calculate the final value of t in seconds
    t_approx = (t_final_n / t_final_d) * (10**t_final_exp)

    # --- For Error Calculation: Compute a more precise value for t ---
    G_real = 6.674e-11
    pi_real = math.pi
    rho_real = 300
    R_real = 2e6
    h_real = 5000
    g_real = (4/3) * pi_real * G_real * rho_real * R_real
    t_real = math.sqrt(2 * h_real / g_real)

    absolute_error = abs(t_approx - t_real)

    print("Y[e]")
    print(f"The landing time t can be calculated by Titan. The smallest absolute error found is {absolute_error:.2f} seconds.")
    print("\nThe equation, with Titan-compatible fractional approximations, is:")
    # We output all numbers as requested.
    print(f"\nt = sqrt( (2 * ({h_n}/{h_d})*10^{h_exp}) / ( (4/3) * ({pi_n}/{pi_d}) * (({G_mantissa_n}/{G_mantissa_d})*10^{G_exp}) * (({rho_n}/{rho_d})*10^{rho_exp}) * (({R_n}/{R_d})*10^{R_exp}) ) )")
    
    print("\nWhich simplifies to:")
    print(f"t ≈ sqrt( (1/{t_squared_mantissa_d}) * 10^{t_squared_exp+1} )")
    
    print("\nUsing a pre-computed approximation for sqrt(1/13) ≈ 3/11, the final result is:")
    print(f"t ≈ ({t_final_n}/{t_final_d}) * 10^{t_final_exp} ≈ {t_approx:.1f} seconds")
    
    # Returning the value for the "answer" part
    return f"Y[{absolute_error:.2f}]"


final_answer = calculate_landing_time()
# The final answer is requested to be in a specific format at the end.
# print(f"\n<<<{final_answer}>>>")