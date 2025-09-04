import math

def check_luminosity_ratio():
    """
    This function checks the correctness of the given answer by recalculating the luminosity ratio
    based on the physical principles described in the problem.
    """

    # --- Define constants and given parameters ---
    # Speed of light in km/s, to match the units of the given velocities.
    c_kms = 299792.458

    # Parameters from the problem statement
    radius_ratio = 1.5  # R1 / R2
    v1_kms = 0.0        # Radial velocity of Star_1 in km/s
    v2_kms = 700.0      # Radial velocity of Star_2 in km/s

    # The answer to check is 'D', which corresponds to a value of ~2.23
    expected_option = 'D'
    expected_value = 2.23

    # --- Recalculate the result based on physics principles ---

    # 1. The Stefan-Boltzmann law gives the luminosity ratio:
    #    L1 / L2 = (R1 / R2)^2 * (T1 / T2)^4
    #    We have (R1 / R2), but we need to find (T1 / T2).

    # 2. Wien's displacement law relates temperature to the rest-frame peak wavelength:
    #    T ∝ 1 / λ_rest
    #    Therefore, T1 / T2 = λ_rest2 / λ_rest1

    # 3. The Doppler shift formula relates the observed wavelength to the rest-frame wavelength:
    #    λ_obs ≈ λ_rest * (1 + v/c)  (for v << c)
    #    So, λ_rest = λ_obs / (1 + v/c)

    # 4. The problem states the observed peak wavelengths are the same: λ_obs1 = λ_obs2.
    #    We can substitute the Doppler formula into the temperature ratio:
    #    T1 / T2 = [λ_obs2 / (1 + v2/c)] / [λ_obs1 / (1 + v1/c)]
    #    Since λ_obs1 = λ_obs2, this simplifies to:
    #    T1 / T2 = (1 + v1/c) / (1 + v2/c)

    temperature_ratio = (1 + v1_kms / c_kms) / (1 + v2_kms / c_kms)

    # 5. Now, substitute the temperature ratio back into the luminosity ratio formula:
    luminosity_ratio = (radius_ratio**2) * (temperature_ratio**4)

    # --- Compare the calculated result with the given answer ---

    # Use a tolerance for floating-point comparison. A 1% tolerance is reasonable for multiple-choice options.
    tolerance = 0.01
    
    if abs(luminosity_ratio - expected_value) < tolerance:
        return "Correct"
    else:
        reason = (
            f"The answer claims the ratio is ~{expected_value} (Option {expected_option}), but the calculation yields a different result.\n"
            f"Constraint check:\n"
            f"1. Radius Ratio (R1/R2): {radius_ratio} -> Used correctly.\n"
            f"2. Velocities (v1, v2): {v1_kms} km/s, {v2_kms} km/s -> Used correctly.\n"
            f"3. Same observed peak wavelength -> This implies T1/T2 = (1+v1/c)/(1+v2/c).\n"
            f"Calculation Steps:\n"
            f" - Temperature Ratio (T1/T2) = (1 + {v1_kms}/{c_kms}) / (1 + {v2_kms}/{c_kms}) = {temperature_ratio:.6f}\n"
            f" - Luminosity Ratio (L1/L2) = (R1/R2)^2 * (T1/T2)^4 = ({radius_ratio})^2 * ({temperature_ratio:.6f})^4 = {luminosity_ratio:.4f}\n"
            f"The calculated ratio {luminosity_ratio:.4f} does not match the expected value of {expected_value} for option {expected_option}."
        )
        return reason

# Execute the check
print(check_luminosity_ratio())