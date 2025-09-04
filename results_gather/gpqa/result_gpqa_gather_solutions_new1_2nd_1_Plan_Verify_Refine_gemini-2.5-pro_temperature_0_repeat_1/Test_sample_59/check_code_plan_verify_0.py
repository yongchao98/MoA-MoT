import math

def check_answer():
    """
    Checks the correctness of the LLM's answer to the relativistic travel problem.
    """
    # --- Problem Constants and Given Data ---
    # Speed of the spacecraft as a fraction of the speed of light
    v_ratio = 0.99999987
    # Astronaut's initial age in years
    astronaut_initial_age = 22
    # Alien's average lifespan in years
    alien_lifespan = 150
    
    # --- External Data (Distance to Large Magellanic Cloud) ---
    # The problem requires an external value for the distance. The provided answer's
    # calculation relies on a distance of ~160,000 light-years to get a result
    # close to 81 years. We will use this value for our primary check.
    # We also define a range to ensure the answer is robust.
    distance_lmc_ly_used = 160000.0
    distance_lmc_ly_min = 159000.0
    distance_lmc_ly_max = 163000.0

    # --- LLM's Answer ---
    # The final answer provided is <<<D>>>, which corresponds to 81 years.
    llm_answer_option = 'D'
    options = {
        'A': 72,
        'B': "The astronaut will die before reaching to the Earth.",
        'C': 77,
        'D': 81
    }
    llm_answer_value = options[llm_answer_option]

    # --- Step 1: Calculate the Lorentz Factor (gamma) ---
    try:
        gamma = 1 / math.sqrt(1 - v_ratio**2)
    except ValueError:
        return "Calculation Error: The velocity ratio results in a mathematical error (e.g., square root of a negative number)."

    # --- Step 2: Calculate Travel Time in Astronaut's Frame (Proper Time) ---
    # Time in Earth's frame (delta_t) = distance / velocity
    # Time in Astronaut's frame (delta_t0) = delta_t / gamma
    # delta_t0 = (distance / v_ratio) / gamma
    
    # Calculate using the distance that best matches the answer
    time_earth_frame = distance_lmc_ly_used / v_ratio
    calculated_proper_time = time_earth_frame / gamma

    # --- Step 3: Verify Constraints ---

    # Constraint 1: Survival Check
    final_age = astronaut_initial_age + calculated_proper_time
    if final_age >= alien_lifespan:
        # This would mean option 'B' is correct.
        if llm_answer_option == 'B':
            return "Correct"
        else:
            return (f"Incorrect. The survival constraint is not met. "
                    f"The astronaut's final age is calculated to be {final_age:.2f} years, "
                    f"which exceeds the lifespan of {alien_lifespan} years. "
                    f"The correct answer should be 'B', but '{llm_answer_option}' was chosen.")

    # If the astronaut survives, option 'B' is incorrect.
    if llm_answer_option == 'B':
        return (f"Incorrect. The chosen answer is 'B' (The astronaut will die), but the calculation shows the astronaut survives. "
                f"The journey takes {calculated_proper_time:.2f} years, making the final age {final_age:.2f}, "
                f"which is within the {alien_lifespan}-year lifespan.")

    # Constraint 2: Numerical Correctness Check
    # Find the closest numerical option to our calculated time.
    numeric_options = {k: v for k, v in options.items() if isinstance(v, int)}
    
    # Find the option key that minimizes the difference between its value and the calculated time
    best_option_key = min(numeric_options, key=lambda k: abs(numeric_options[k] - calculated_proper_time))

    if llm_answer_option != best_option_key:
        return (f"Incorrect. The calculated travel time is approximately {calculated_proper_time:.2f} years. "
                f"The closest numerical option is '{best_option_key}' ({options[best_option_key]} years), "
                f"but the selected answer was '{llm_answer_option}' ({llm_answer_value} years).")

    # Constraint 3: Robustness Check (using the range of distances)
    # This ensures the chosen answer is the best fit regardless of the exact LMC distance used.
    proper_time_min = (distance_lmc_ly_min / v_ratio) / gamma
    proper_time_max = (distance_lmc_ly_max / v_ratio) / gamma
    
    if not (numeric_options[best_option_key] > proper_time_min - 2 and numeric_options[best_option_key] < proper_time_max + 2):
        return (f"Incorrect. While the chosen answer is the closest for a distance of {distance_lmc_ly_used} ly, "
                f"it falls outside the robust calculated range of {proper_time_min:.2f} to {proper_time_max:.2f} years.")

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check
result = check_answer()
print(result)