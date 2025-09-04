import math

def check_answer():
    """
    This function checks the correctness of the provided answer by recalculating the result from the problem description.
    """
    # --- Problem Parameters ---
    # Effective temperature of the star in Kelvin
    T_eff = 6000.0
    # Temperature difference of the spots in Kelvin
    delta_T = 1000.0
    # Filling factor of the spots on one hemisphere
    f = 0.20

    # --- Candidate Answer ---
    # The final answer provided is 'A', which corresponds to the value ~0.32
    # from the question's option list.
    expected_value = 0.32

    # --- Step 1: Calculate the amplitude of brightness variation from starspots ---
    # The flux (F) is proportional to the fourth power of temperature (T), F ∝ T⁴.
    # The amplitude is the relative drop in flux when the spotted hemisphere is visible.
    # Amplitude = f * (1 - (T_spot / T_eff)⁴)

    # Calculate the spot temperature
    T_spot = T_eff - delta_T

    # Calculate the amplitude
    try:
        amplitude_spots = f * (1 - (T_spot / T_eff)**4)
    except Exception as e:
        return f"An error occurred during the amplitude calculation: {e}"

    # --- Step 2: Relate amplitude to an exoplanet transit ---
    # The transit depth (amplitude) is the square of the planet-to-star radius ratio.
    # Transit Depth = (R_pl / R_star)²
    # To have the same signal, we equate the two amplitudes:
    # (R_pl / R_star)² = amplitude_spots

    # --- Step 3: Solve for the relative radius (R_pl / R_star) ---
    # R_pl / R_star = sqrt(amplitude_spots)
    try:
        calculated_ratio = math.sqrt(amplitude_spots)
    except ValueError:
        return f"Calculation error: Cannot take the square root of a negative amplitude ({amplitude_spots})."
    except Exception as e:
        return f"An error occurred during the square root calculation: {e}"

    # --- Step 4: Compare the calculated result with the expected answer ---
    # We use math.isclose() to account for potential rounding in the options.
    # A relative tolerance of 2% is reasonable for "approximately equal to".
    if math.isclose(calculated_ratio, expected_value, rel_tol=0.02):
        return "Correct"
    else:
        reason = (
            f"Incorrect. The provided answer corresponds to a radius ratio of {expected_value}, "
            f"but the calculated ratio is approximately {calculated_ratio:.4f}.\n\n"
            "Here is the step-by-step calculation:\n"
            f"1. Spot Temperature (T_spot) = {T_eff} K - {delta_T} K = {T_spot} K.\n"
            "2. The amplitude of the spot-induced variation is calculated as:\n"
            f"   Amplitude = f * (1 - (T_spot / T_eff)⁴)\n"
            f"   Amplitude = {f} * (1 - ({T_spot} / {T_eff})⁴) ≈ {amplitude_spots:.5f}\n"
            "3. The equivalent exoplanet transit depth is (R_pl / R_star)², so we solve for the ratio:\n"
            f"   R_pl / R_star = sqrt(Amplitude)\n"
            f"   R_pl / R_star = sqrt({amplitude_spots:.5f}) ≈ {calculated_ratio:.5f}\n"
            f"The calculated value {calculated_ratio:.5f} is not consistent with the answer's value of {expected_value}."
        )
        return reason

# Execute the check and print the result
print(check_answer())