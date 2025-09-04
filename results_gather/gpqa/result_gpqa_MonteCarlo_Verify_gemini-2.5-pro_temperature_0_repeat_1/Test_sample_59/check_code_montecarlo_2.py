import math

def check_relativity_problem():
    """
    This function checks the correctness of the provided answer to the special relativity problem.
    It calculates the proper time for the astronaut and compares it with the given options and constraints.
    """
    # --- Define constants and given values ---
    # Speed as a fraction of the speed of light (v/c)
    v_over_c = 0.99999987
    # Distance to the Large Magellanic Cloud in light-years
    distance_light_years = 163000
    # Astronaut's initial age
    astronaut_initial_age = 22
    # Alien average lifespan
    alien_lifespan = 150

    # The provided answer key and the options from the question
    provided_answer_key = 'D'
    options = {
        'A': 77,
        'B': 72,
        'C': 'The astronaut will die before reaching to the Earth.',
        'D': 81
    }

    # --- Step 1: Calculate time from Earth's perspective (Δt) ---
    # For a distance in light-years and speed as a fraction of c, the time in years is d/v.
    # Δt = distance_light_years / v_over_c
    time_earth = distance_light_years / v_over_c

    # --- Step 2: Calculate the Lorentz factor (γ) ---
    # γ = 1 / sqrt(1 - (v/c)^2)
    try:
        gamma = 1 / math.sqrt(1 - v_over_c**2)
    except ValueError:
        return "Calculation Error: Speed (v) cannot be greater than or equal to the speed of light (c)."

    # --- Step 3: Calculate the time experienced by the astronaut (proper time, Δt₀) ---
    # Δt₀ = Δt / γ
    time_astronaut = time_earth / gamma

    # --- Step 4: Check constraints and correctness of the answer ---
    remaining_lifespan = alien_lifespan - astronaut_initial_age
    survives = time_astronaut < remaining_lifespan

    # Check if the survival status matches the provided answer
    if provided_answer_key == 'C':
        if survives:
            return (f"Incorrect. The calculated travel time for the astronaut is {time_astronaut:.2f} years, "
                    f"which is less than their remaining lifespan of {remaining_lifespan} years. "
                    f"The astronaut survives, but the answer is C (the astronaut dies).")
        else:
            # The astronaut dies and the answer is C, which is consistent.
            return "Correct"
    else: # The provided answer is A, B, or D
        if not survives:
            return (f"Incorrect. The calculated travel time for the astronaut is {time_astronaut:.2f} years, "
                    f"which is more than their remaining lifespan of {remaining_lifespan} years. "
                    f"The astronaut would die, so the correct answer should be C.")

    # If the astronaut survives and the answer is not C, find the closest numerical option.
    numerical_options = {k: v for k, v in options.items() if isinstance(v, (int, float))}
    
    # Find the key of the option with the minimum absolute difference from the calculated time
    closest_option_key = min(
        numerical_options.keys(),
        key=lambda k: abs(numerical_options[k] - time_astronaut)
    )

    # Compare the closest option with the provided answer
    if closest_option_key == provided_answer_key:
        return "Correct"
    else:
        return (f"Incorrect. The calculated travel time for the astronaut is {time_astronaut:.2f} years. "
                f"The closest option is {closest_option_key} ({options[closest_option_key]} years), "
                f"but the provided answer was {provided_answer_key} ({options[provided_answer_key]} years).")

# Run the check
result = check_relativity_problem()
print(result)