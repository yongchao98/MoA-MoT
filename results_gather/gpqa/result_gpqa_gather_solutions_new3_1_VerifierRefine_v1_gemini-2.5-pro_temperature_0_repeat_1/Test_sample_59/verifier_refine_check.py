import math

def check_correctness():
    """
    Checks the correctness of the LLM's answer to the special relativity problem.
    """
    # --- Problem Constraints & Given Data ---
    speed_ratio = 0.99999987  # v/c
    lifespan = 150
    initial_age = 22

    # --- Options from the question ---
    # The final answer provided is <<<A>>>, which corresponds to 81 years.
    options = {'A': 81, 'C': 72, 'D': 77}
    llm_answer_value = options['A']

    # --- External Data: Plausible distances to the LMC in light-years ---
    # This is the main source of variation in the final calculation.
    plausible_distances = {
        "Lower Estimate": 159000,
        "Commonly Rounded": 160000,
        "Upper Estimate (Nature 2013)": 163000
    }

    # --- Physics Calculation ---
    # Calculate the Lorentz factor (gamma)
    try:
        gamma = 1 / math.sqrt(1 - speed_ratio**2)
    except ValueError:
        return "Calculation Error: The speed ratio must be less than 1."

    # --- Verification Loop ---
    for source, distance_ly in plausible_distances.items():
        # Calculate time in Earth's frame (t_earth)
        t_earth = distance_ly / speed_ratio

        # Calculate time for the astronaut (proper time, t_astronaut)
        t_astronaut = t_earth / gamma

        # 1. Check the survival constraint
        remaining_lifespan = lifespan - initial_age
        if t_astronaut >= remaining_lifespan:
            return (f"Incorrect. The survival constraint is not met for a distance of {distance_ly} ly. "
                    f"The journey takes {t_astronaut:.2f} years, but the astronaut only has "
                    f"{remaining_lifespan} years of life left.")

        # 2. Check if the calculated time is closest to the selected answer (81 years)
        min_diff = float('inf')
        closest_option = None
        for key, value in options.items():
            diff = abs(t_astronaut - value)
            if diff < min_diff:
                min_diff = diff
                closest_option = value
        
        if closest_option != llm_answer_value:
            return (f"Incorrect. For a distance of {distance_ly} ly ({source}), the calculated time is "
                    f"{t_astronaut:.2f} years. This is closest to {closest_option} years, "
                    f"not the selected answer of {llm_answer_value} years.")

    # If all checks pass for all plausible distances, the answer is robust and correct.
    return "Correct"

# Run the check and print the result
result = check_correctness()
print(result)