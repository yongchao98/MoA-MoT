import math

def check_transit_probability():
    """
    This function checks the correctness of the answer to the exoplanet transit probability question.

    The problem states:
    - P2 = 3 * P1 (Orbital period of Planet_2 is three times that of Planet_1)
    - M_s1 = 2 * M_s2 (Mass of host star 1 is twice that of host star 2)
    - R_s1 = R_s2 (Radii of host stars are the same)

    The transit probability (p_tr) is proportional to R_s / a.
    The ratio of probabilities is p_tr1 / p_tr2 = (R_s1 / a1) / (R_s2 / a2).
    Since R_s1 = R_s2, this simplifies to p_tr1 / p_tr2 = a2 / a1.

    From Kepler's Third Law, a^3 is proportional to M_s * P^2.
    So, a is proportional to (M_s * P^2)^(1/3).
    The ratio a2 / a1 = [ (M_s2 * P2^2) / (M_s1 * P1^2) ]^(1/3).
    """

    # Let's use relative values to perform the calculation.
    # Let Planet_2's system have base values.
    P2 = 3.0
    M_s2 = 1.0

    # Calculate Planet_1's system properties based on the given relationships.
    P1 = P2 / 3.0
    M_s1 = 2.0 * M_s2

    # Calculate the ratio of the semi-major axes, a2 / a1.
    # The proportionality constant cancels out in the ratio.
    # a2 / a1 = [ (M_s2 * P2^2) / (M_s1 * P1^2) ]^(1/3)
    ratio_a2_over_a1 = ((M_s2 * P2**2) / (M_s1 * P1**2))**(1/3)

    # The ratio of transit probabilities p_tr1 / p_tr2 is equal to a2 / a1.
    prob_ratio_p1_over_p2 = ratio_a2_over_a1

    # The provided answer is 'B', which states:
    # "Planet_1 is preferred due to its ~1.65 times higher probability to transit."
    
    # Constraint 1: Planet_1 must be preferred.
    # This is true if its probability is higher, i.e., the ratio > 1.
    if prob_ratio_p1_over_p2 <= 1:
        return f"Incorrect. The answer states Planet_1 is preferred, but the calculated probability ratio (P1/P2) is {prob_ratio_p1_over_p2:.4f}, which is not greater than 1."

    # Constraint 2: The ratio must be approximately 1.65.
    expected_ratio = 1.65
    if not math.isclose(prob_ratio_p1_over_p2, expected_ratio, rel_tol=0.01):
        return f"Incorrect. The answer states the probability ratio is ~{expected_ratio}, but the calculated ratio is {prob_ratio_p1_over_p2:.4f}."

    # If all constraints are satisfied, the answer is correct.
    return "Correct"

# Run the check
result = check_transit_probability()
print(result)