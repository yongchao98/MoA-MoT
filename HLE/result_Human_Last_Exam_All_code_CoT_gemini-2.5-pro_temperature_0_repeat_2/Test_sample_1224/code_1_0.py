import math

def solve_pandora_landing_time():
    """
    Calculates the landing time on Pandora using Titan's 4-bit architecture constraints
    and determines the calculation error.
    """

    # --- 1. Define constants and parameters ---
    h = 5000.0  # m

    # Simplified model parameters for Titan
    R_titan = 2 * 10**6  # m
    d_titan = 300.0      # kg/m^3

    # Approximations for Titan computation
    # pi is approximated as 3/1
    pi_num, pi_den = 3, 1
    # G is approximated as 2/3 * 10^-10
    G_num, G_den = 2, 3
    G_exp = -10

    # --- 2. Perform the Titan-feasible calculation ---
    # The formula is t^2 = (3/2) * (h/R) * (1/d) * (1/pi) * (1/G)
    # We calculate the mantissa and exponent of t^2 separately.

    # Term 1: 3/2
    t1_num, t1_den = 3, 2
    t1_exp = 0

    # Term 2: h/R
    t2_num, t2_den = 5, 2 # from 5000 / (2 * 10^6) = 2.5 * 10^-3 = 5/2 * 10^-3
    t2_exp = -3

    # Term 3: 1/d
    t3_num, t3_den = 1, 3 # from 1 / 300 = 1/3 * 10^-2
    t3_exp = -2

    # Term 4: 1/pi
    t4_num, t4_den = pi_den, pi_num
    t4_exp = 0

    # Term 5: 1/G
    t5_num, t5_den = G_den, G_num
    t5_exp = -G_exp

    # Multiply mantissas step-by-step, simplifying at each stage
    # This demonstrates feasibility as all intermediate numerators/denominators are <= 15
    # (3/2) * (5/2) = 15/4
    # (15/4) * (1/3) = 15/12 = 5/4
    # (5/4) * (1/3) = 5/12
    # (5/12) * (3/2) = 15/24 = 5/8
    t_squared_mant_num = 5
    t_squared_mant_den = 8

    # Sum exponents
    t_squared_exp = t1_exp + t2_exp + t3_exp + t4_exp + t5_exp

    # Combine mantissa and exponent for t^2
    # t^2 = 5/8 * 10^5 = 0.625 * 10^5 = 6.25 * 10^4
    # This can also be written as 25/4 * 10^4, which is a perfect square
    t_squared_val = (t_squared_mant_num / t_squared_mant_den) * (10**t_squared_exp)

    # Calculate t
    # t = sqrt(25/4 * 10^4) = 5/2 * 10^2 = 250
    t_titan = math.sqrt(t_squared_val)

    # --- 3. Calculate the 'true' value for error comparison ---
    G_true = 6.67430e-11
    pi_true = math.pi
    R_true = 2e6 # Using the same simplified model for a fair comparison
    d_true = 300.0

    g_true = (4/3) * pi_true * G_true * R_true * d_true
    t_true = math.sqrt(2 * h / g_true)

    # --- 4. Calculate the absolute error ---
    error = abs(t_titan - t_true)

    # --- 5. Print the results ---
    print("Titan Calculation Steps:")
    print("Formula: t^2 = (3/2) * (h/R) * (1/d) * (1/pi) * (1/G)")
    print(f"Approximations: pi = {pi_num}/{pi_den}, G = {G_num}/{G_den} x 10^{G_exp}")
    print("\nFinal Equation on Titan:")
    # The final equation shows all the numbers used in the calculation
    final_equation = (
        f"t = sqrt( ( {t1_num}/{t1_den} ) * "
        f"( ({t2_num}*10^{3}) / ({t2_den}*10^{6}) ) * "
        f"( 1 / ({t3_den}*10^{2}) ) * "
        f"( {t4_den}/{t4_num} ) * "
        f"( {t5_num}/{t5_den}*10^{t5_exp} ) )"
    )
    print(final_equation)
    print(f"\nTitan Calculated Result: t = {t_titan:.2f} s")
    print(f"More Precise Result: t = {t_true:.2f} s")
    print(f"Smallest Absolute Error: e = {error:.2f}")
    
    # The final answer format
    print(f"\n<<<Y[{error:.2f}]>>>")

solve_pandora_landing_time()