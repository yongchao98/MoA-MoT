import math

def check_answer():
    """
    Checks the correctness of the provided answer for the physics problem.

    The problem involves two parts:
    1. Calculating the fraction of maximum power radiated at a specific angle.
    2. Determining the dependence of radiated power on the wavelength.

    The physical model is an oscillating electric dipole.
    """

    # The LLM's chosen answer is C.
    # Let's define the options from the question.
    options = {
        'A': {'fraction': 1/4, 'lambda_exponent': -3},
        'B': {'fraction': 3/4, 'lambda_exponent': -6},
        'C': {'fraction': 1/4, 'lambda_exponent': -4},
        'D': {'fraction': 1/2, 'lambda_exponent': -4}
    }
    chosen_option_key = 'C'
    chosen_option = options[chosen_option_key]

    # --- Part 1: Verify the fraction of power at 30 degrees ---
    # For an oscillating electric dipole, the radiated power per unit solid angle (dP/dΩ)
    # is proportional to sin^2(theta), where theta is the angle with the oscillation axis.
    # dP/dΩ = K * sin^2(theta)

    # The maximum power 'A' occurs at theta = 90 degrees (pi/2 radians), where sin^2(90) = 1.
    # So, A is proportional to 1.
    max_power_factor = math.sin(math.radians(90))**2

    # We need to find the power at theta = 30 degrees.
    angle_deg = 30
    angle_rad = math.radians(angle_deg)
    power_factor_at_30 = math.sin(angle_rad)**2

    # The fraction of the maximum power is the ratio of the power factors.
    calculated_fraction = power_factor_at_30 / max_power_factor

    # --- Part 2: Verify the wavelength dependence ---
    # The power (P) radiated by an electric dipole is proportional to the fourth power
    # of the angular frequency (omega). P ∝ ω^4.
    # This is a fundamental result from Larmor's formula and its relativistic generalization.
    # Angular frequency (ω) is related to wavelength (λ) by ω = 2*π*c/λ, so ω ∝ 1/λ.
    # Therefore, P ∝ (1/λ)^4 = λ^(-4).
    # The exponent of lambda is -4.
    calculated_lambda_exponent = -4

    # --- Part 3: Compare calculated results with the chosen option ---
    expected_fraction = chosen_option['fraction']
    expected_lambda_exponent = chosen_option['lambda_exponent']

    # Check if both parts of the answer are correct
    fraction_is_correct = math.isclose(calculated_fraction, expected_fraction, rel_tol=1e-9)
    exponent_is_correct = (calculated_lambda_exponent == expected_lambda_exponent)

    if fraction_is_correct and exponent_is_correct:
        return "Correct"
    else:
        error_messages = []
        if not fraction_is_correct:
            error_messages.append(
                f"The fraction of maximum power is incorrect. "
                f"Based on the sin^2(theta) dependence, the fraction at 30 degrees should be {calculated_fraction:.4f}, "
                f"but the chosen answer states it is {expected_fraction}."
            )
        if not exponent_is_correct:
            error_messages.append(
                f"The wavelength dependence is incorrect. "
                f"For electric dipole radiation, power is proportional to lambda^({calculated_lambda_exponent}), "
                f"but the chosen answer states it is lambda^({expected_lambda_exponent})."
            )
        return "Incorrect. " + " ".join(error_messages)

# Execute the check and print the result
result = check_answer()
print(result)