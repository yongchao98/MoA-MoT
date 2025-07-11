import numpy as np

def solve_stellar_measurement():
    """
    This script checks the consistency of stellar measurements using Planck's Law
    and identifies the most likely source of error if an inconsistency is found.
    """

    # Constants
    h = 6.62607015e-34  # Planck constant in J*s
    c = 299792458      # Speed of light in m/s
    k_B = 1.380649e-23   # Boltzmann constant in J/K

    # Given measurements
    L_nm = 400.0         # Wavelength in nanometers
    L_m = L_nm * 1e-9    # Wavelength in meters
    B_measured = 1.2e15  # Measured spectral radiance in W/(m^2*sr*m)
    T_measured = 9000.0    # Measured surface temperature in Kelvin

    # First, we check if the measured radiance is physically possible for the given temperature.
    # The maximum radiance for a given temperature occurs at the peak wavelength,
    # determined by Wien's Displacement Law.
    lambda_peak = (2.898e-3) / T_measured # Wien's constant in m*K
    
    # Calculate the maximum possible radiance at T_measured
    def planck_radiance(wav, temp):
        term1 = (2 * h * c**2) / (wav**5)
        exponent = (h * c) / (wav * k_B * temp)
        term2 = np.exp(exponent) - 1
        return term1 / term2

    B_max_possible = planck_radiance(lambda_peak, T_measured)

    if B_measured > B_max_possible:
        # The measured radiance is higher than the theoretical maximum for T=9000K.
        # This confirms an error in either T or B. We assume T is the error.
        # We rearrange Planck's law to solve for the correct temperature.
        # T = hc / (λ * k_B * ln( (2 * h * c^2) / (λ^5 * B) + 1 ) )
        arg_ln = (2 * h * c**2) / (L_m**5 * B_measured) + 1
        T_correct = (h * c) / (L_m * k_B * np.log(arg_ln))
        T_correct_rounded = int(round(T_correct))

        print("An inconsistency was found in the measurements.")
        print(f"The measured radiance {B_measured:.2e} is greater than the maximum possible radiance ({B_max_possible:.2e}) for a star at {T_measured} K.")
        print("Assuming the Temperature (T) measurement is the error, the corrected value is calculated.\n")
        
        # Print the final equation with the corrected temperature
        print("Planck's Law with the corrected temperature and original measurements:")
        print(f"{B_measured:.1e} = (2 * {h:.4e} * {c:.4e}^2) / "
              f"({L_m:.1e}^5 * (exp(({h:.4e} * {c:.4e}) / "
              f"({L_m:.1e} * {k_B:.4e} * {T_correct_rounded})) - 1))")
        
        # Verify that the right side of the equation now equals the measured radiance
        B_recalculated = planck_radiance(L_m, T_correct_rounded)
        print(f"\nVerification: The right side of the equation evaluates to approximately {B_recalculated:.2e}.")
        print(f"This is consistent with the measured spectral radiance.")

        # Final answer
        answer = f"T{T_correct_rounded}"
        print(f"\nConclusion: The quantity most likely in error is T, and its corrected value is {T_correct_rounded} K.")

    else:
        # This block would execute if the measurements were consistent, but they are not.
        answer = "0"
        print("The measurements appear to be consistent.")
        print(f"The result is: {answer}")


solve_stellar_measurement()