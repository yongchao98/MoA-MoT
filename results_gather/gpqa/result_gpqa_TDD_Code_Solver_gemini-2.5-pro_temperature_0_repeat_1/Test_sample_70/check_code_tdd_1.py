import math

def check_exoplanet_temperature_ratio():
    """
    This function verifies the solution to the exoplanet temperature ratio problem.

    The problem asks for the ratio of equilibrium temperatures between Planet_4 and Planet_2 (T_eq,4 / T_eq,2).

    The key physical relationships are:
    1. Kepler's Third Law: For planets orbiting the same star, the semi-major axis 'a' and orbital period 'P' are related by a^3 ∝ P^2. This implies a ∝ P^(2/3).
    2. Equilibrium Temperature: For planets with the same albedo orbiting the same star, the equilibrium temperature T_eq is proportional to 1/sqrt(a).

    Combining these, we get:
    T_eq,4 / T_eq,2 = (a_2 / a_4)^(1/2)
                   = [(P_2 / P_4)^(2/3)]^(1/2)
                   = (P_2 / P_4)^(1/3)
    """

    # Given data from the question
    period_ratios = {'Planet_1': 1, 'Planet_2': 2, 'Planet_3': 2.5, 'Planet_4': 3.5, 'Planet_5': 5}
    P2 = period_ratios['Planet_2']
    P4 = period_ratios['Planet_4']

    # The options provided in the question
    options = {
        "A": 0.75,
        "B": 0.57,
        "C": 0.83,
        "D": 0.69
    }
    
    # The answer given by the LLM
    llm_answer_label = "C"

    # --- Step 1: Calculate the correct value based on the physics ---
    try:
        # The formula is (P_2 / P_4)^(1/3)
        calculated_ratio = (P2 / P4)**(1/3)
    except (ZeroDivisionError, KeyError) as e:
        return f"Error in calculation setup: {e}"

    # --- Step 2: Check if the LLM's answer is the closest option ---
    
    # Find the option that is numerically closest to our calculated value
    closest_option_label = min(options, key=lambda k: abs(options[k] - calculated_ratio))

    # --- Step 3: Compare and return the result ---
    if llm_answer_label == closest_option_label:
        # Further check if the value is reasonably close to the option's value
        if math.isclose(calculated_ratio, options[llm_answer_label], rel_tol=0.02): # Allow 2% tolerance for rounding
             return "Correct"
        else:
            return (f"Incorrect. The calculated ratio is {calculated_ratio:.4f}. "
                    f"While '{llm_answer_label}' is the closest option, the value {options[llm_answer_label]} "
                    f"is not a good approximation of the calculated result.")
    else:
        return (f"Incorrect. The calculated temperature ratio is (P2/P4)^(1/3) = ({P2}/{P4})^(1/3) ≈ {calculated_ratio:.4f}. "
                f"This value is closest to option '{closest_option_label}' ({options[closest_option_label]}), "
                f"not the provided answer '{llm_answer_label}' ({options[llm_answer_label]}).")

# Run the check
result = check_exoplanet_temperature_ratio()
print(result)