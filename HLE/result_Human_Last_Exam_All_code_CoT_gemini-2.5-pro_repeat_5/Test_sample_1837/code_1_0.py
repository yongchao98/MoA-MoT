import math

def solve_star_measurement():
    """
    Analyzes the consistency of stellar measurements using Planck's Law.
    """
    # Physical constants
    h = 6.62607015e-34  # Planck's constant (J*s)
    c = 2.99792458e8    # Speed of light (m/s)
    k = 1.380649e-23    # Boltzmann constant (J/K)

    # Measured values
    lambda_val = 400e-9  # Wavelength in meters (400 nm)
    b_meas = 1.2e15      # Spectral radiance in W/(m^2*sr*m)
    t_meas = 9000        # Temperature in Kelvin

    # --- Consistency Check ---
    # Calculate theoretical spectral radiance B(T, L)
    try:
        exponent_val = (h * c) / (lambda_val * k * t_meas)
        denominator = math.exp(exponent_val) - 1
        numerator = 2 * h * c**2
        b_calc = numerator / (lambda_val**5 * denominator)
    except (OverflowError, ZeroDivisionError):
        b_calc = float('inf')

    # --- Analysis ---
    # The calculated radiance (b_calc, ~2.18e14) is much lower than the
    # measured radiance (b_meas, 1.2e15). The measurements are not consistent.
    # Answering "0" is not correct. We must find the likely error.

    # Scenario 1: Is the wavelength L wrong?
    # The maximum radiance for T=9000K occurs at lambda_peak ≈ 322 nm.
    # The radiance at this peak is ~2.39e14, which is still less than b_meas.
    # Therefore, no wavelength at T=9000K can produce the measured radiance.
    # The error is not in L.

    # Scenario 2: Is the temperature T wrong?
    # Solve Planck's law for T, assuming B and L are correct.
    # T = hc / (L*k * ln( (2hc^2 / (L^5 * B)) + 1 ))
    try:
        term_in_log = (2 * h * c**2) / (lambda_val**5 * b_meas) + 1
        if term_in_log <= 0:
            # This case implies b_meas is too high, which we already established
            # would be impossible if T was correct. But here we assume T is wrong.
            # A negative log argument means no real solution for T.
            # Let's re-check math.
            # B = C1 / (L^5 * (exp(C2/LT)-1)) -> exp(C2/LT)-1 = C1/(B*L^5)
            # exp(C2/LT) = C1/(B*L^5) + 1. Since C1, B, L are positive, the arg is > 1.
            # So a real solution must exist.
            pass

        log_term = math.log(term_in_log)
        t_calc = (h * c) / (lambda_val * k * log_term)
        t_calc_rounded = round(t_calc)

    except (ValueError, ZeroDivisionError):
        t_calc_rounded = None

    # Conclusion: A corrected temperature of ~15180 K is highly plausible for
    # a DB-class white dwarf, whereas the initially measured 9000 K is low.
    # This makes temperature the most likely source of error.

    print("The inconsistency suggests a measurement error. Based on physical plausibility, the temperature (T) is the most likely incorrect quantity.")
    print("Assuming the spectral radiance and wavelength are correct, we can calculate the expected temperature using this equation:")
    
    # Print the equation with all the numbers
    equation_str = (
        f"T = ({h:.5e} * {c:.5e}) / "
        f"({lambda_val:.1e} * {k:.5e} * "
        f"ln((2 * {h:.5e} * ({c:.5e})^2) / (({lambda_val:.1e})^5 * {b_meas:.1e}) + 1))"
    )
    print(equation_str)

    answer = f"T{t_calc_rounded}"
    print(f"\nThe calculated correct temperature is {t_calc_rounded} K.")
    print(f'<<<{answer}>>>')

solve_star_measurement()