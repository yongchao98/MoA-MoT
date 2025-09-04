import math

def check_correctness():
    """
    Checks the correctness of the answer to the exoplanet transit probability question.

    The final answer to check is 'B', which states:
    "Planet_1 is preferred due to its ~1.65 times higher probability to transit."
    """

    # 1. Define the relationships from the problem statement as ratios.
    # P_2 = 3 * P_1  => P_ratio_P2_over_P1 = 3
    P_ratio_P2_over_P1 = 3.0
    # M_s1 = 2 * M_s2 => M_ratio_Ms1_over_Ms2 = 2
    M_ratio_Ms1_over_Ms2 = 2.0
    # R_s1 = R_s2 => R_ratio_Rs1_over_Rs2 = 1
    R_ratio_Rs1_over_Rs2 = 1.0

    # 2. Calculate the ratio of probabilities p1 / p2 based on physical laws.
    # The derivation is as follows:
    # p_transit ∝ R_s / a
    # a ∝ (P^2 * M_s)^(1/3)
    # p1 / p2 = (R_s1 / a1) / (R_s2 / a2) = (R_s1 / R_s2) * (a2 / a1)
    # a2 / a1 = [ (P_2^2 * M_s2) / (P_1^2 * M_s1) ]^(1/3)
    # a2 / a1 = [ (P_2/P_1)^2 * (M_s2/M_s1) ]^(1/3)
    # a2 / a1 = [ P_ratio_P2_over_P1^2 * (1 / M_ratio_Ms1_over_Ms2) ]^(1/3)

    # Calculate the ratio of semi-major axes a2/a1
    try:
        a2_over_a1 = ( (P_ratio_P2_over_P1**2) * (1 / M_ratio_Ms1_over_Ms2) )**(1/3)
    except Exception as e:
        return f"An error occurred during calculation: {e}"

    # The probability ratio p1/p2 is R_ratio * a2/a1
    prob_ratio_p1_over_p2 = R_ratio_Rs1_over_Rs2 * a2_over_a1

    # 3. Define the expected outcome from the given answer 'B'.
    expected_preferred_planet = "Planet_1"
    expected_ratio = 1.65

    # 4. Verify the calculated results against the expected outcome.
    # Check which planet is preferred
    if prob_ratio_p1_over_p2 > 1:
        calculated_preferred_planet = "Planet_1"
    elif prob_ratio_p1_over_p2 < 1:
        calculated_preferred_planet = "Planet_2"
    else:
        calculated_preferred_planet = "Neither (equal probability)"

    # Check if the preferred planet matches the answer
    if calculated_preferred_planet != expected_preferred_planet:
        return (f"Incorrect: The preferred planet is wrong. "
                f"The answer states {expected_preferred_planet} is preferred, but the calculation "
                f"shows {calculated_preferred_planet} is preferred. The calculated probability "
                f"ratio p1/p2 is {prob_ratio_p1_over_p2:.4f}.")

    # Check if the ratio factor matches the answer
    # We use math.isclose for robust floating-point comparison, with a 1.5% relative tolerance
    # to account for the "~" (approximately) in the option.
    if not math.isclose(prob_ratio_p1_over_p2, expected_ratio, rel_tol=0.015):
        return (f"Incorrect: The probability ratio is wrong. "
                f"The answer states the factor is ~{expected_ratio}, but the calculated "
                f"factor is {prob_ratio_p1_over_p2:.4f}.")

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check and print the result.
result = check_correctness()
print(result)