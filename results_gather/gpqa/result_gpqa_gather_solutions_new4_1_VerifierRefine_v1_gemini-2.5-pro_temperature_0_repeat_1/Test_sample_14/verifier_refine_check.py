import math

def check_answer():
    """
    Checks the correctness of the answer to the exoplanet transit probability question.
    """
    # --- Define relationships from the problem statement ---
    # Let's use relative values for simplicity.
    # Let P1 = 1, then P2 = 3.
    # Let Ms2 = 1, then Ms1 = 2.
    # Let Rs1 = 1, Rs2 = 1.
    
    # P_ratio = P2 / P1
    P_ratio = 3.0
    # M_ratio = Ms1 / Ms2
    M_ratio = 2.0
    # R_ratio = Rs1 / Rs2
    R_ratio = 1.0

    # --- Derivation of the probability ratio ---
    # Transit probability (p_tr) is proportional to Rs / a
    # p_tr1 / p_tr2 = (Rs1 / a1) / (Rs2 / a2) = (Rs1 / Rs2) * (a2 / a1)
    
    # From Kepler's Third Law: P^2 is proportional to a^3 / Ms
    # So, a is proportional to (P^2 * Ms)^(1/3)
    
    # The ratio of semi-major axes a2 / a1 is:
    # a2 / a1 = [ (P2^2 * Ms2) / (P1^2 * Ms1) ]^(1/3)
    # a2 / a1 = [ (P2/P1)^2 / (Ms1/Ms2) ]^(1/3)
    
    try:
        # Calculate the ratio of semi-major axes (a2 / a1)
        a2_over_a1 = (P_ratio**2 / M_ratio)**(1/3)
        
        # The ratio of probabilities p_tr1 / p_tr2 is (Rs1/Rs2) * (a2/a1)
        # Since Rs1/Rs2 = 1, the probability ratio is equal to a2/a1
        prob_ratio_1_to_2 = R_ratio * a2_over_a1

    except Exception as e:
        return f"An error occurred during calculation: {e}"

    # --- Check against the provided answer ---
    # The final answer is 'A', which states:
    # "Planet_1 is preferred due to its ~1.65 times higher probability to transit."
    
    expected_ratio = 1.65
    expected_preferred_planet = "Planet_1"
    
    # Check which planet is preferred based on the calculation
    if prob_ratio_1_to_2 > 1:
        calculated_preferred_planet = "Planet_1"
    elif prob_ratio_1_to_2 < 1:
        calculated_preferred_planet = "Planet_2"
    else:
        calculated_preferred_planet = "Neither"
        
    # Check if the calculated preferred planet matches the answer
    if calculated_preferred_planet != expected_preferred_planet:
        return f"Incorrect: The answer states {expected_preferred_planet} is preferred, but the calculation shows {calculated_preferred_planet} is preferred because the probability ratio p1/p2 is {prob_ratio_1_to_2:.2f}."

    # Check if the calculated ratio matches the answer's ratio
    # We use math.isclose for floating point comparison, with a tolerance.
    # A tolerance of 0.01 is reasonable for "~1.65".
    if not math.isclose(prob_ratio_1_to_2, expected_ratio, rel_tol=0.01):
        return f"Incorrect: The answer states the probability ratio is ~{expected_ratio}, but the calculation gives a ratio of {prob_ratio_1_to_2:.4f}."

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check
result = check_answer()
print(result)