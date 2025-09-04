import math

def check_answer():
    """
    Checks the correctness of the final answer for the relativistic travel problem.
    """
    # --- Problem Parameters ---
    v_over_c = 0.99999987  # Speed of the spacecraft as a fraction of c
    astronaut_initial_age = 22  # years
    alien_lifespan = 150  # years

    # The distance to the Large Magellanic Cloud (LMC) is not given.
    # We use a range of accepted astronomical values.
    # The value of 160,000 light-years is commonly used and aligns best with the options.
    distance_ly_best_fit = 160000.0  # light-years
    distance_ly_range = (159000.0, 163000.0) # A plausible range for the distance

    # --- Question Options ---
    # Note: The final answer uses the original question's option lettering.
    options = {
        "A": 77,
        "B": "The astronaut will die before reaching to the Earth.",
        "C": 72,
        "D": 81
    }
    final_answer_key = "D"
    final_answer_value = options[final_answer_key]

    # --- Step 1: Calculate the Lorentz Factor (gamma) ---
    try:
        gamma = 1 / math.sqrt(1 - v_over_c**2)
    except ValueError:
        return "Incorrect: Calculation error. The term inside the square root is negative, which is physically impossible for v < c."

    # --- Step 2: Calculate Travel Time for a range of distances ---
    # Time in Earth's frame (delta_t)
    delta_t_min = distance_ly_range[0] / v_over_c
    delta_t_max = distance_ly_range[1] / v_over_c
    
    # Time in Astronaut's frame (delta_t0)
    delta_t0_min = delta_t_min / gamma
    delta_t0_max = delta_t_max / gamma

    # --- Step 3: Verify the Survival Constraint ---
    # The astronaut's age upon arrival will be between:
    final_age_min = astronaut_initial_age + delta_t0_min
    final_age_max = astronaut_initial_age + delta_t0_max

    if final_age_min >= alien_lifespan:
        # This means the astronaut would die regardless of the exact distance.
        if final_answer_key != "B":
            return f"Incorrect: The astronaut's final age would be at least {final_age_min:.2f} years, which exceeds their lifespan of {alien_lifespan} years. The answer should be B."
        else:
            # The answer correctly states the astronaut dies.
            pass # Continue to check if this is the most likely outcome
    
    if final_answer_key == "B" and final_age_max < alien_lifespan:
        return f"Incorrect: The provided answer is B (The astronaut will die), but the calculated final age is between {final_age_min:.2f} and {final_age_max:.2f} years, which is within the {alien_lifespan}-year lifespan. The astronaut survives."

    # --- Step 4: Find the best numerical option ---
    # Use the "best fit" distance for a precise comparison
    delta_t_best_fit = distance_ly_best_fit / v_over_c
    delta_t0_best_fit = delta_t_best_fit / gamma

    numerical_options = {k: v for k, v in options.items() if isinstance(v, int)}
    if not numerical_options:
        return "Incorrect: No numerical options to compare against."

    # Find the option closest to our calculated best-fit time
    closest_option_key = min(
        numerical_options.keys(),
        key=lambda k: abs(numerical_options[k] - delta_t0_best_fit)
    )

    # --- Step 5: Final Verification ---
    if final_answer_key != closest_option_key:
        return (f"Incorrect: The provided answer is {final_answer_key} ({final_answer_value} years), but the calculation shows the most plausible answer is {closest_option_key} ({options[closest_option_key]} years). "
                f"The calculated travel time is approximately {delta_t0_best_fit:.2f} years, which is closest to {options[closest_option_key]}.")

    # All checks passed
    return "Correct"

# Run the check
result = check_answer()
print(result)