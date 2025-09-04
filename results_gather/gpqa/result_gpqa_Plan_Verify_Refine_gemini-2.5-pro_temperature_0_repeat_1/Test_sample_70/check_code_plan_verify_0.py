import math

def check_answer():
    """
    This function verifies the solution to the exoplanet temperature ratio problem.
    
    The key physical relationships are:
    1. Equilibrium Temperature (T_eq) is proportional to 1 / sqrt(a), where 'a' is the semi-major axis.
       This comes from balancing the absorbed stellar flux with the radiated thermal energy.
       Absorbed ~ 1/a^2, Radiated ~ T_eq^4. So T_eq^4 ~ 1/a^2, which means T_eq ~ 1/sqrt(a).
    2. Kepler's Third Law states that the square of the orbital period (P) is proportional to the cube of the semi-major axis (a).
       P^2 ~ a^3, which means a ~ P^(2/3).
    
    Combining these, we get:
    T_eq ~ 1 / sqrt(P^(2/3)) => T_eq ~ 1 / P^(1/3).
    
    Therefore, the ratio of equilibrium temperatures between Planet_4 and Planet_2 is:
    T_eq4 / T_eq2 = (P2 / P4)^(1/3).
    """
    
    # Constraint: Orbital period ratios
    period_ratios = {
        "Planet_1": 1,
        "Planet_2": 2,
        "Planet_3": 2.5,
        "Planet_4": 3.5,
        "Planet_5": 5
    }
    
    # Extract the relevant periods for the ratio T_eq4 / T_eq2
    try:
        P2 = period_ratios["Planet_2"]
        P4 = period_ratios["Planet_4"]
    except KeyError as e:
        return f"Constraint Error: The required planet {e} is not in the provided period ratios."

    # Calculate the theoretical temperature ratio
    # T_eq4 / T_eq2 = (P2 / P4)^(1/3)
    calculated_ratio = (P2 / P4)**(1/3)
    
    # The LLM's answer is 'A', which corresponds to a value of ~0.83
    llm_answer_option = 'A'
    options = {
        'A': 0.83,
        'B': 0.75,
        'C': 0.69,
        'D': 0.57
    }
    
    # Find the closest option to the calculated value
    # This confirms if the LLM selected the correct option based on the calculation
    closest_option = min(options, key=lambda k: abs(options[k] - calculated_ratio))
    
    if closest_option != llm_answer_option:
        return (f"The answer is incorrect. The calculated ratio is {calculated_ratio:.4f}, "
                f"which is closest to option {closest_option} ({options[closest_option]}), "
                f"not option {llm_answer_option}.")

    # Check if the value of the chosen option is a reasonable approximation
    # We use a relative tolerance of 2% to account for the "~" (approximately) sign.
    if not math.isclose(calculated_ratio, options[llm_answer_option], rel_tol=0.02):
        return (f"The answer is technically correct in choosing the closest option, but the approximation is poor. "
                f"Calculated value is {calculated_ratio:.4f}, while option {llm_answer_option} is {options[llm_answer_option]}.")

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check
result = check_answer()
print(result)