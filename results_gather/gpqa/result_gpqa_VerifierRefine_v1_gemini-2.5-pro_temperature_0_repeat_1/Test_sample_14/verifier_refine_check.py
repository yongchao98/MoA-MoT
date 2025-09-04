import math

def check_transit_probability_answer():
    """
    This function checks the correctness of the provided answer to the exoplanet transit probability problem.
    It calculates the ratio of transit probabilities for Planet_1 and Planet_2 based on the given physical relationships
    and compares it to the selected answer (Option D).

    The key formulas used are:
    1. Transit Probability: P ≈ R_star / a
       - R_star is the stellar radius.
       - a is the orbital semi-major axis.
    2. Kepler's Third Law: T^2 / a^3 ∝ 1 / M_star
       - T is the orbital period.
       - M_star is the stellar mass.
    """
    try:
        # --- Step 1: Define the relationships from the problem statement as ratios ---
        # The orbital period of Planet_1 is three times shorter than that of Planet_2.
        # T1 = T2 / 3  => T1/T2 = 1/3
        T1_div_T2 = 1.0 / 3.0
        
        # The star hosting Planet_1 has a mass that is twice that of the host star of Planet_2.
        # M_s1 = 2 * M_s2 => M_s1/M_s2 = 2
        M_s1_div_M_s2 = 2.0
        
        # Both host stars have the same radii.
        # R_s1 = R_s2 => R_s1/R_s2 = 1
        R_s1_div_R_s2 = 1.0

        # --- Step 2: Calculate the ratio of semi-major axes (a2 / a1) ---
        # From Kepler's Third Law, we can derive the relationship for the semi-major axis:
        # a^3 ∝ M_s * T^2  =>  a ∝ (M_s * T^2)^(1/3)
        #
        # We can then find the ratio of a1 to a2:
        # a1/a2 = [ (M_s1/M_s2) * (T1/T2)^2 ]^(1/3)
        a1_div_a2 = (M_s1_div_M_s2 * (T1_div_T2)**2)**(1.0/3.0)
        
        # For the probability ratio, we need the inverse, a2/a1.
        a2_div_a1 = 1.0 / a1_div_a2

        # --- Step 3: Calculate the ratio of transit probabilities (P1 / P2) ---
        # The transit probability P is proportional to R_s / a.
        # The ratio of probabilities P1 / P2 is:
        # P1/P2 = (R_s1 / a1) / (R_s2 / a2) = (R_s1/R_s2) * (a2/a1)
        P1_div_P2_calculated = R_s1_div_R_s2 * a2_div_a1

        # --- Step 4: Analyze the provided answer (Option D) ---
        # The answer is D: "Planet_1 is preferred due to its ~1.65 times higher probability to transit."
        expected_preferred_planet = "Planet_1"
        expected_ratio = 1.65
        # A relative tolerance of 2% is reasonable for the "~" (approximately) sign.
        relative_tolerance = 0.02 

        # --- Step 5: Perform the checks to verify the answer ---
        
        # Check 1: Verify the preferred planet.
        # If P1/P2 > 1, Planet_1 has a higher probability and is preferred.
        if P1_div_P2_calculated > 1:
            calculated_preferred_planet = "Planet_1"
        elif P1_div_P2_calculated < 1:
            calculated_preferred_planet = "Planet_2"
        else:
            calculated_preferred_planet = "Neither"

        if calculated_preferred_planet != expected_preferred_planet:
            return (f"Incorrect: The choice of the preferred planet is wrong. "
                    f"The answer states {expected_preferred_planet} is preferred, but the calculation shows "
                    f"{calculated_preferred_planet} is preferred because the probability ratio P1/P2 is "
                    f"{P1_div_P2_calculated:.4f}, which is {'greater' if P1_div_P2_calculated > 1 else 'less'} than 1.")

        # Check 2: Verify the numerical ratio.
        if not math.isclose(P1_div_P2_calculated, expected_ratio, rel_tol=relative_tolerance):
            return (f"Incorrect: The probability ratio is wrong. "
                    f"The answer states the ratio is ~{expected_ratio}, but the calculated ratio is "
                    f"{P1_div_P2_calculated:.4f}, which is not within a {relative_tolerance*100}% tolerance of the expected value.")

        # If all checks pass, the answer is correct.
        return "Correct"

    except Exception as e:
        return f"An error occurred during the checking process: {e}"

# Execute the check and print the result.
result = check_transit_probability_answer()
print(result)