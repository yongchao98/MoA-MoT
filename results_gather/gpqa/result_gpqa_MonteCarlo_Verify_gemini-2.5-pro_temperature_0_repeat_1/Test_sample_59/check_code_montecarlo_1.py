import math

def check_relativity_answer():
    """
    This function checks the correctness of the answer to the special relativity problem.
    """
    # --- Problem Constants ---
    # Speed of the spacecraft as a fraction of the speed of light (c)
    v_fraction = 0.99999987
    # Astronaut's age at the start of the journey in solar years
    astronaut_age = 22
    # Average lifespan of the alien species in solar years
    alien_lifespan = 150
    # The most commonly cited distance to the Large Magellanic Cloud in light-years.
    # This is a necessary piece of external information to solve the problem.
    distance_ly = 163000.0
    
    # The answer provided by the LLM
    llm_answer = 'D'

    # --- Calculations ---
    
    # 1. Calculate the remaining lifespan of the astronaut
    remaining_lifespan = alien_lifespan - astronaut_age
    if remaining_lifespan <= 0:
        return "Invalid input: Astronaut's age is greater than or equal to the lifespan."

    # 2. Calculate the Lorentz factor (gamma)
    try:
        gamma = 1 / math.sqrt(1 - v_fraction**2)
    except ValueError:
        return "Calculation Error: The fraction of the speed of light must be less than 1."

    # 3. Calculate the travel time as measured by a stationary observer (e.g., on Earth)
    # Since distance is in light-years and speed is a fraction of c, the time will be in years.
    # time_earth = distance_ly / (v_fraction * c), where c = 1 light-year/year
    time_earth_frame = distance_ly / v_fraction

    # 4. Calculate the time experienced by the astronaut (proper time) due to time dilation
    time_astronaut_frame = time_earth_frame / gamma

    # --- Verification ---

    # Constraint 1: Check if the astronaut survives the journey
    if time_astronaut_frame > remaining_lifespan:
        # The astronaut would die. The correct answer should be C.
        if llm_answer == 'C':
            return "Correct"
        else:
            return (f"Incorrect. The calculated travel time for the astronaut is {time_astronaut_frame:.2f} years. "
                    f"This exceeds the astronaut's remaining lifespan of {remaining_lifespan} years. "
                    f"Therefore, the correct answer is C, but the provided answer was {llm_answer}.")

    # If the astronaut survives, find the closest numerical option.
    options = {'A': 77, 'B': 72, 'D': 81}
    
    # Calculate the absolute difference between the calculated time and each option
    differences = {key: abs(value - time_astronaut_frame) for key, value in options.items()}
    
    # Find the option with the smallest difference
    closest_option = min(differences, key=differences.get)

    # Final check: Does the LLM's answer match the closest calculated option?
    if llm_answer == closest_option:
        return "Correct"
    else:
        return (f"Incorrect. The calculated travel time for the astronaut is {time_astronaut_frame:.2f} years. "
                f"This value is closest to option {closest_option} ({options[closest_option]} years). "
                f"The provided answer was {llm_answer}.")

# Execute the check and print the result
result = check_relativity_answer()
print(result)