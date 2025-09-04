import math

def check_answer():
    """
    This function verifies the solution to the particle physics problem.

    The problem involves two scenarios of particle decay, and we need to find the
    Lorentz factor for the second scenario.

    The fraction of surviving particles (N/N0) is given by the formula:
    N/N0 = exp(-t_lab / (gamma * tau))
    where:
    - t_lab is the time of flight in the lab frame (R/c).
    - gamma is the Lorentz factor.
    - tau is the proper mean lifetime of the particle.

    The term (R/c)/tau is a constant for a given particle and detector.
    Let K = (R/c)/tau. The formula becomes N/N0 = exp(-K / gamma).
    By taking the natural log, we get: ln(N/N0) = -K / gamma.
    This can be rearranged to solve for K or gamma.
    """

    # --- Scenario 1 ---
    # Given Lorentz factor
    gamma_1 = 20
    # Given fraction of surviving particles
    fraction_1 = 1/3

    # --- Scenario 2 ---
    # Desired fraction of surviving particles
    fraction_2 = 2/3
    # The proposed answer from the LLM is C, which corresponds to a Lorentz factor of 54.
    proposed_answer_value = 54
    proposed_answer_option = 'C'

    # --- Calculation ---
    # From the formula ln(N/N0) = -K / gamma, we can express K for both scenarios:
    # K = -gamma_1 * ln(fraction_1)
    # K = -gamma_2 * ln(fraction_2)
    #
    # Setting them equal gives:
    # -gamma_1 * ln(fraction_1) = -gamma_2 * ln(fraction_2)
    #
    # Solving for gamma_2:
    # gamma_2 = gamma_1 * ln(fraction_1) / ln(fraction_2)
    #
    # Using logarithm properties ln(1/x) = -ln(x) and ln(x/y) = -ln(y/x):
    # gamma_2 = gamma_1 * (-ln(3)) / (-ln(3/2))
    # gamma_2 = gamma_1 * ln(3) / ln(3/2)

    try:
        # Calculate the required Lorentz factor for the second scenario
        calculated_gamma_2 = gamma_1 * math.log(3) / math.log(1.5)
    except Exception as e:
        return f"An error occurred during calculation: {e}"

    # --- Verification ---
    # The options are integers, so we check if the proposed answer is the closest integer
    # to the calculated value. A small tolerance is used to account for potential
    # floating-point inaccuracies, but checking for the nearest integer is more robust.
    if round(calculated_gamma_2) == proposed_answer_value:
        return "Correct"
    else:
        reason = (
            f"The provided answer is {proposed_answer_option} ({proposed_answer_value}), but the calculation yields a different result.\n"
            f"The derivation is correct: γ₂ = γ₁ * ln(3) / ln(3/2).\n"
            f"Using γ₁ = 20, the calculated value for γ₂ is approximately {calculated_gamma_2:.2f}.\n"
            f"The closest integer to the calculated value is {round(calculated_gamma_2)}, which does not match the proposed answer of {proposed_answer_value}."
        )
        return reason

# Run the check
result = check_answer()
print(result)