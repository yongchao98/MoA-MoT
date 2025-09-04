import math

def check_answer():
    """
    Checks the correctness of the LLM's answer for the physics problem.
    """
    # The LLM's answer is B, which corresponds to a coefficient of 0.506 for lambda/a.
    llm_answer_coefficient = 0.506

    # --- Theoretical Calculation ---
    # The problem simplifies to Fraunhofer diffraction from a circular aperture of radius 'a'.
    # The angular distance between the first two minima is given by:
    # delta_theta = ((z2 - z1) / (2 * pi)) * (lambda / a)
    # where z1 and z2 are the first two non-zero roots of the Bessel function J1(x).

    # Values for the first two non-zero roots of J1(x).
    # Using high-precision values.
    z1 = 3.83170597
    z2 = 7.01558667

    # Calculate the theoretical coefficient for (lambda / a).
    calculated_coefficient = (z2 - z1) / (2 * math.pi)

    # --- Verification ---
    # Compare the calculated coefficient with the one from the LLM's answer.
    # The options are given to three decimal places, so a relative tolerance of 1% is appropriate.
    if math.isclose(calculated_coefficient, llm_answer_coefficient, rel_tol=1e-2):
        return "Correct"
    else:
        reason = (
            f"The answer is incorrect. The provided answer corresponds to a coefficient of {llm_answer_coefficient}, "
            f"but the calculated coefficient is approximately {calculated_coefficient:.4f}.\n"
            f"Calculation details:\n"
            f"  - The problem describes diffraction from a circular aperture of radius 'a'.\n"
            f"  - The angular distance between the first two minima is Δθ = ((z₂ - z₁) / (2π)) * (λ/a).\n"
            f"  - Using z₁ ≈ {z1:.4f} and z₂ ≈ {z2:.4f}, the coefficient is ({z2:.4f} - {z1:.4f}) / (2π) ≈ {calculated_coefficient:.4f}."
        )
        return reason

# Run the check
result = check_answer()
print(result)