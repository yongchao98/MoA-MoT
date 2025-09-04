import math

def check_relativity_answer():
    """
    This function checks the correctness of the LLM's answer to the relativity problem.
    It verifies the physics calculations and the logical conclusion based on those calculations.
    """
    # --- Problem Constraints and Given Values ---
    v_over_c = 0.99999987  # Speed of the spacecraft as a fraction of the speed of light
    alien_lifespan = 150  # years
    astronaut_start_age = 22  # years

    # --- Assumption made in the LLM's answer ---
    # The distance to the Large Magellanic Cloud is assumed to be 160,000 light-years,
    # a common approximation that leads to one of the multiple-choice options.
    distance_ly = 160000  # light-years

    # --- LLM's chosen answer ---
    llm_answer_choice = 'C'
    llm_answer_value = 81  # years

    # --- Step 1: Calculate the Lorentz factor (gamma) ---
    # gamma = 1 / sqrt(1 - (v/c)^2)
    try:
        gamma = 1 / math.sqrt(1 - v_over_c**2)
    except ValueError:
        return "Error: Calculation for gamma resulted in a math domain error. This can happen if v >= c."

    # --- Step 2: Calculate the time in Earth's reference frame (delta_t) ---
    # For a distance in light-years, time in years = distance_ly / (v/c)
    time_earth_frame = distance_ly / v_over_c

    # --- Step 3: Calculate the time for the astronaut (proper time, delta_t0) ---
    # According to time dilation, delta_t0 = delta_t / gamma
    time_astronaut_frame = time_earth_frame / gamma

    # --- Step 4: Verify the survival constraint ---
    # The astronaut must complete the journey before their lifespan ends.
    remaining_lifespan = alien_lifespan - astronaut_start_age
    if time_astronaut_frame > remaining_lifespan:
        return (f"Incorrect. The astronaut would die before reaching Earth. "
                f"The journey takes {time_astronaut_frame:.2f} years from the astronaut's perspective, "
                f"but they only have {remaining_lifespan} years of life remaining. "
                f"The correct answer should be D, not C.")

    # --- Step 5: Verify if the chosen answer is the closest approximation ---
    # The question asks for an approximate time. We check if the calculated time is
    # closest to the value of the chosen option C.
    options = {'A': 72, 'B': 77, 'C': 81}
    
    # Calculate the absolute difference between the calculated time and each option's value
    differences = {opt: abs(time_astronaut_frame - val) for opt, val in options.items()}
    
    # Find the option with the smallest difference
    closest_option_key = min(differences, key=differences.get)

    if closest_option_key == llm_answer_choice:
        # The calculation confirms that 81 years is the closest answer and the astronaut survives.
        # The reasoning and conclusion of the LLM are sound.
        # print(f"Calculation successful: Astronaut's travel time is ~{time_astronaut_frame:.2f} years. Closest option is C.")
        return "Correct"
    else:
        return (f"Incorrect. The chosen answer C is not the closest option. "
                f"Based on a distance of {distance_ly} light-years, the calculated travel time is "
                f"{time_astronaut_frame:.2f} years. The closest option to this value is "
                f"{closest_option_key} ({options[closest_option_key]} years), not C.")

# Run the check
result = check_relativity_answer()
print(result)