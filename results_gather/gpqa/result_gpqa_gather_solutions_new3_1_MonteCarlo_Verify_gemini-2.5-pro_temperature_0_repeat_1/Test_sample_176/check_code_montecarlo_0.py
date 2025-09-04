import math

def check_luminosity_ratio():
    """
    Calculates the luminosity ratio based on the problem's parameters
    and checks if it matches the provided answer 'C' (~2.23).
    """
    # --- Given parameters ---
    # Radius ratio of Star_1 to Star_2
    radius_ratio = 1.5
    # Radial velocity of Star_2 in km/s
    v2_kms = 700.0
    # Speed of light in km/s (approximation used in the problem's solutions)
    c_kms = 300000.0

    # --- Physics Calculation ---
    # The luminosity ratio L1/L2 = (R1/R2)^2 * (T1/T2)^4

    # 1. Calculate the radius term
    radius_term = radius_ratio ** 2
    if not math.isclose(radius_term, 2.25):
        return "Constraint check failed: The radius ratio squared (1.5^2) should be 2.25."

    # 2. Calculate the temperature ratio term
    # From Doppler effect and Wien's law, T1/T2 = 1 / (1 + v2/c)
    v_over_c = v2_kms / c_kms
    temperature_ratio = 1 / (1 + v_over_c)
    temperature_term = temperature_ratio ** 4

    # 3. Calculate the final luminosity ratio
    calculated_ratio = radius_term * temperature_term

    # --- Verification ---
    # The provided answer is 'C', which corresponds to ~2.23
    expected_value = 2.23
    
    # Check if the calculated value rounds to the expected value
    if round(calculated_ratio, 2) != expected_value:
        return (f"Incorrect: The calculated luminosity ratio is {calculated_ratio:.4f}, "
                f"which rounds to {round(calculated_ratio, 2)}. This does not match the "
                f"expected value of {expected_value} for answer 'C'.")

    # Check if the reasoning in the provided answer is sound.
    # The reasoning correctly identifies all physical principles and combines them.
    # It correctly notes that the mass information is extraneous.
    # It correctly interprets "appeared brightest" as the observed wavelength,
    # necessitating the use of the Doppler effect.
    # The final calculation L₁ / L₂ ≈ 2.229, which rounds to 2.23, is correct.
    
    return "Correct"

# Run the check
result = check_luminosity_ratio()
print(result)