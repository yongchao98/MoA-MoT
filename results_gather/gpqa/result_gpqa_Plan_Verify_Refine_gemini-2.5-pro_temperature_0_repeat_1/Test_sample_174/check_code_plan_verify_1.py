import numpy as np

def check_dipole_radiation_answer():
    """
    This function checks the correctness of the provided answer for the
    oscillating spheroid (electric dipole) radiation problem.

    It verifies two physical constraints:
    1. The fraction of maximum power radiated at θ = 30°.
    2. The dependence of radiated power on the wavelength λ.
    """

    # --- Constraint 1: Power Fraction at θ = 30° ---
    # The radiated power per unit solid angle is dP/dΩ ∝ sin²(θ).
    # The maximum power, A, occurs when sin²(θ) is maximum, which is 1 (at θ = 90°).
    # We need to find the power at θ = 30° as a fraction of A.
    # Fraction = (Power at 30°) / (Maximum Power) = sin²(30°) / sin²(90°) = sin²(30°).
    
    theta_deg = 30
    theta_rad = np.deg2rad(theta_deg)
    
    # Calculate the theoretically correct fraction
    correct_fraction = np.sin(theta_rad)**2

    # --- Constraint 2: Wavelength Dependence ---
    # The radiated power is proportional to the fourth power of the angular frequency (ω⁴).
    # The relationship between angular frequency ω and wavelength λ is ω = 2πc/λ.
    # Therefore, the power is proportional to (1/λ)⁴ = λ⁻⁴.
    correct_exponent = -4

    # --- Extract values from the chosen answer (Option C) ---
    # The question options are:
    # A) 1/4, λ^(-3)
    # B) 1/2, λ^(-4)
    # C) 1/4, λ^(-4)
    # D) 3/4, λ^(-6)
    # The LLM's answer is C.
    answer_fraction = 1/4
    answer_exponent = -4

    # --- Verification ---
    # Check if the fraction from the answer matches the calculated correct fraction.
    # We use np.isclose for safe floating-point comparison.
    is_fraction_correct = np.isclose(correct_fraction, answer_fraction)
    
    # Check if the exponent from the answer matches the correct exponent.
    is_exponent_correct = (correct_exponent == answer_exponent)

    if is_fraction_correct and is_exponent_correct:
        return "Correct"
    else:
        error_messages = []
        if not is_fraction_correct:
            error_messages.append(
                f"Constraint on power fraction is not satisfied. "
                f"The radiated power is proportional to sin²(θ). At θ=30°, the fraction of maximum power should be sin²(30°) = (0.5)² = {correct_fraction:.2f}. "
                f"The chosen answer's fraction is {answer_fraction}."
            )
        if not is_exponent_correct:
            error_messages.append(
                f"Constraint on wavelength dependence is not satisfied. "
                f"The radiated power is proportional to ω⁴, which means it is proportional to λ^({correct_exponent}). "
                f"The chosen answer's wavelength dependence is λ^({answer_exponent})."
            )
        return "\n".join(error_messages)

# Run the check
result = check_dipole_radiation_answer()
print(result)