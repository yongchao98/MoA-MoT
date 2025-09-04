import math

def check_answer():
    """
    Checks the correctness of the given answer for the cobalt-thiocyanate complex problem.
    """
    # --- Given parameters from the question ---
    ligand_conc = 0.1  # [SCN-] in M
    betas = [9, 40, 63, 16]  # [β1, β2, β3, β4]
    n_of_interest = 2  # For the dithiocyanato complex [Co(SCN)2]

    # The answer to check is C, which is 16.9%
    expected_percentage = 16.9

    # --- Calculation ---

    # Calculate the denominator (Φ), which is the sum of relative concentrations of all species.
    # Φ = 1 (for free Co^2+) + β1*[L]^1 + β2*[L]^2 + β3*[L]^3 + β4*[L]^4
    phi = 1.0
    for i, beta in enumerate(betas):
        power = i + 1
        term = beta * (ligand_conc ** power)
        phi += term

    # Calculate the numerator for the species of interest, [Co(SCN)2]
    # Numerator = β2 * [L]^2
    # Note: betas list is 0-indexed, so β2 is at index 1.
    numerator = betas[n_of_interest - 1] * (ligand_conc ** n_of_interest)

    # Calculate the fraction (α_n) and convert to percentage
    if phi == 0:
        calculated_percentage = 0.0
    else:
        fraction = numerator / phi
        calculated_percentage = fraction * 100

    # --- Verification ---
    # Check if the calculated percentage is close to the expected answer (16.9%)
    # We use a tolerance to account for floating-point arithmetic and rounding in the option.
    if math.isclose(calculated_percentage, expected_percentage, rel_tol=1e-2):
        return "Correct"
    else:
        return (f"Incorrect. The calculated percentage is {calculated_percentage:.2f}%. "
                f"The provided answer C corresponds to {expected_percentage}%, which is not consistent with the calculation. "
                f"The calculated value {calculated_percentage:.2f}% is the correct percentage for the dithiocyanato cobalt(II) complex.")

# Run the check
result = check_answer()
print(result)