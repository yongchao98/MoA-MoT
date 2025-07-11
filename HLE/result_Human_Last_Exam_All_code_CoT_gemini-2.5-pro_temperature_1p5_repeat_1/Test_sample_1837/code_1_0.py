import numpy as np

def solve_pandora_case():
    """
    Checks the consistency of stellar measurements using Planck's Law and identifies the most likely error.
    """
    # 1. Define physical constants and measured values
    h = 6.62607015e-34  # Planck's constant in J*s
    c = 299792458       # Speed of light in m/s
    k_B = 1.380649e-23   # Boltzmann constant in J/K

    # Measured values from the user
    lambda_nm = 400
    lambda_m = lambda_nm * 1e-9  # Wavelength in meters
    T_measured = 9000             # Temperature in Kelvin
    B_measured = 1.2e15           # Spectral radiance in W/(m^2 * sr * m)

    # Helper function to calculate spectral radiance using Planck's Law
    def planck_radiance(lam_m, temp_k):
        """Calculates Planck's black-body spectral radiance."""
        exponent = (h * c) / (lam_m * k_B * temp_k)
        # Avoid overflow for large exponents, though not expected here
        if exponent > 709:
            return 0.0
        return (2 * h * c**2) / (lam_m**5 * (np.exp(exponent) - 1))

    # 2. Check if the measured radiance is physically possible for the given temperature.
    # The peak of the Planck curve for a given temperature can be found using the root of d(B)/d(lambda)=0.
    # The solution for x = hc/(lambda*k*T) is approximately 4.965114.
    wien_const_x = 4.965114231744276
    lambda_peak_m = (h * c) / (wien_const_x * k_B * T_measured)
    B_max_at_T_measured = planck_radiance(lambda_peak_m, T_measured)

    print(f"Checking consistency for measurements: T={T_measured} K, L={lambda_nm} nm, B={B_measured:.1e} W/(m^2*sr*m)")
    print(f"The theoretical maximum spectral radiance for a body at {T_measured} K is {B_max_at_T_measured:.2e} W/(m^2*sr*m).")
    
    # 3. If B_measured > B_max, an error in T or B is confirmed.
    if B_measured > B_max_at_T_measured:
        print("\nThe measured radiance is higher than the theoretical maximum for the given temperature.")
        print("This implies a measurement error in either temperature (T) or spectral radiance (B).")
        print("Given that 9000 K is unusually cool for a DB-class white dwarf, the temperature reading is the most likely error.")
        
        # 4. Calculate the corrected temperature T_corrected, assuming L and B are correct.
        print("\nCalculating the corrected temperature by solving Planck's Law for T:")
        print("Formula: T = (h*c) / (L * k_B * ln((2*h*c^2 / (L^5 * B)) + 1))")
        
        # Breaking down the calculation for display
        numerator = h * c
        term_in_log_num = 2 * h * c**2
        term_in_log_den = lambda_m**5 * B_measured
        log_argument = (term_in_log_num / term_in_log_den) + 1
        denominator = lambda_m * k_B * np.log(log_argument)
        
        T_corrected = numerator / denominator
        T_corrected_rounded = int(round(T_corrected))

        # 5. Output each number in the final equation
        print("\nPlugging in the numbers for the corrected temperature calculation:")
        print(f"T_corrected = ({numerator:.4e}) / ({lambda_m:.1e} * {k_B:.4e} * ln(({term_in_log_num:.4e} / (({lambda_m:.1e})^5 * {B_measured:.1e})) + 1))")
        print(f"T_corrected = {numerator:.4e} / ({denominator:.4e})")
        print(f"Calculated T_corrected = {T_corrected:.1f} K")
        print(f"Rounding to the nearest unit gives a corrected temperature of {T_corrected_rounded} K.")

        result = f"T{T_corrected_rounded}"
    else:
        # This case is for when the measurements are physically possible, but might still be inconsistent.
        # Here we check if the provided values satisfy Planck's Law within a tolerance.
        B_expected = planck_radiance(lambda_m, T_measured)
        if np.isclose(B_measured, B_expected, rtol=0.05): # using a 5% relative tolerance
             print("\nThe measurements are physically consistent.")
             result = "0"
        else:
            # If not consistent, we assume error in B as T is more robustly determined than flux.
            print("\nMeasurements are inconsistent. Assuming error in spectral radiance B.")
            B_expected_rounded = int(round(B_expected))
            result = f"B{B_expected_rounded}"

    print(f"<<<{result}>>>")

solve_pandora_case()