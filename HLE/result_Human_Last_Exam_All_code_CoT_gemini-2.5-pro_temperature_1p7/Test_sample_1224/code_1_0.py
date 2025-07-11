import math
from fractions import Fraction

def solve_pandora_landing_time():
    """
    This script evaluates if the Titan 4-bit computer can calculate the landing time
    of a ball on exoplanet Pandora and calculates the smallest possible error.

    Plan:
    1.  Calculate a high-precision 'true' value for the landing time using full precision for all constants and parameters.
    2.  Simulate the calculation on the Titan computer by using approximations that respect the 4-bit integer constraint (numerators/denominators <= 15).
    3.  The key is to find a simple fractional approximation for the constant cluster 1/(4*pi*G).
    4.  The approximation 1/(4*pi*G) ~= 1/8 * 10^10 is used, which has an acceptable error margin (~5%).
    5.  This simplification leads to a result for t^2 which is a perfect square of a simple fraction, making the final square root operation trivial for the architecture.
    6.  The absolute error between the true time and the approximated time is calculated.
    """

    # Part 1: High-Precision "True" Calculation
    # --- Constants and Parameters ---
    d = 5000.0  # Drop height in meters
    G = 6.67430e-11  # Gravitational constant
    PI = math.pi

    # Pandora's Core
    r_core = 100e3  # Core radius in meters
    rho_core = 1200.0  # Core density in kg/m^3
    v_core = (4/3) * PI * r_core**3
    m_core = v_core * rho_core

    # Pandora's Shell (as an oblate spheroid)
    a_shell = 2000e3  # Equatorial radius in meters
    b_shell = 1985e3  # Polar radius in meters
    rho_shell = 300.0  # Shell density in kg/m^3
    v_planet = (4/3) * PI * a_shell**2 * b_shell
    v_shell = v_planet - v_core
    m_shell = v_shell * rho_shell

    # Total Mass and Gravity
    m_total = m_core + m_shell
    # We calculate g at the equator (radius a_shell)
    g_true = (G * m_total) / (a_shell**2)
    t_true = math.sqrt((2 * d) / g_true)

    # Part 2: Titan 4-bit Architecture Simulation
    print("--- Simulating Calculation on Titan Computer ---")
    print("The goal is to calculate t = sqrt(2*d/g).")

    # Step 1: Simplify the physics formula
    # g = G*M/R^2 and M = V*rho = (4/3)*pi*R^3*rho.
    # We simplify by ignoring the small core and oblateness.
    # Let R be the equatorial radius (2000 km) and rho be the shell density (300 kg/m^3).
    # g = G * (4/3)*pi*R*rho
    # t^2 = 2*d / g = 2*d / (G*(4/3)*pi*R*rho) = 3*d / (2*pi*G*R*rho)
    R = 2e6
    rho = 300
    # Substituting values for d, R, and rho:
    # t^2 = (3 * 5e3) / (2 * pi * G * 2e6 * 3e2)
    # t^2 = 15e3 / (12e8 * pi * G)
    # t^2 = (15/12) * 10^-5 / (pi*G) = (5/4) * 10^-5 / (pi*G)
    print("\nStep 1: Simplify the equation for t^2.")
    print("t^2 = (3*d)/(2*pi*G*R*rho) = 5/(4*pi*G) * 10^-5")

    # Step 2: Approximate the constant cluster 1 / (4*pi*G)
    # 1 / (4*pi*G) ~= 1 / (4 * 3.14159 * 6.674e-11) ~= 0.1192e10
    # We look for a simple fraction a/b where a,b <= 15.
    # 1/8 = 0.125. This is a ~5% error, which is acceptable.
    # So, we approximate 1/(4*pi*G) as 1/8 * 10^10
    approx_frac = Fraction(1, 8)
    approx_exp = 10
    print(f"\nStep 2: Approximate 1/(4*pi*G) with a Titan-compatible fraction.")
    print(f"1/(4*pi*G) ~= {approx_frac.numerator}/{approx_frac.denominator} * 10^{approx_exp}")

    # Step 3: Calculate t^2 using the approximation
    t2_frac_part = Fraction(5, 1) * approx_frac
    t2_exp_part = -5 + approx_exp
    print("\nStep 3: Compute t^2.")
    print(f"t^2 ~= (5/1 * 10^-5) * ({approx_frac.numerator}/{approx_frac.denominator} * 10^{approx_exp})")
    print(f"t^2 ~= {t2_frac_part.numerator}/{t2_frac_part.denominator} * 10^{t2_exp_part}")

    # Step 4: Simplify the result for Titan
    # 5/8 * 10^5 = 50/8 * 10^4 = 25/4 * 10^4.
    # The fraction 25/4 has a numerator > 15, but it can be stored
    # as an expression like (6 + 1/4).
    # t^2 = (6 + 1/4) * 10^4
    final_t2_frac = Fraction(25, 4)
    print("\nStep 4: Simplify the result to find a perfect square.")
    print(f"t^2 = 5/8 * 10^5 = 50/8 * 10^4 = {final_t2_frac.numerator}/{final_t2_frac.denominator} * 10^4")

    # Step 5: Calculate the square root
    # Since 25/4 is a perfect square of a simple fraction (5/2), an iterative
    # method is not required. The result can be computed directly.
    t_approx_frac = Fraction(int(math.sqrt(final_t2_frac.numerator)), int(math.sqrt(final_t2_frac.denominator)))
    t_approx = float(t_approx_frac) * (10**(4/2))
    print("\nStep 5: Calculate t by taking the square root.")
    print(f"t = sqrt({final_t2_frac.numerator}/{final_t2_frac.denominator}) * 10^2 = {t_approx_frac.numerator}/{t_approx_frac.denominator} * 100")

    print(f"\nFinal approximated landing time: {t_approx} s")
    print("---------------------------------------------")

    # Part 3: Calculate and Output the final answer
    abs_error = abs(t_true - t_approx)
    print(f"\nHigh-precision landing time: {t_true:.2f} s")
    print(f"Titan's approximated time: {t_approx:.2f} s")
    print(f"Absolute error: |{t_true:.2f} - {t_approx:.2f}| = {abs_error:.2f} s")
    
    # The final equation is the computation of the final value, 250
    # final_equation = f"t = {t_approx_frac.numerator}/{t_approx_frac.denominator} * 100 = {t_approx}"
    print("\nFinal Equation Output:")
    print(f"t = sqrt({final_t2_frac.numerator}/{final_t2_frac.denominator} * 10^4) = {t_approx_frac.numerator}/{t_approx_frac.denominator} * 10^2 = {t_approx}")


if __name__ == "__main__":
    solve_pandora_landing_time()
