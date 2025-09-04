import math

def check_correctness():
    """
    Checks the correctness of the final answer to the special relativity problem.
    The final answer provided is 'C', which corresponds to 81 years in the original question.
    """
    
    # --- Problem Constraints and Given Data ---
    # The final answer to check, corresponding to option C.
    final_answer_value = 81  # years

    # Speed of the spacecraft as a fraction of the speed of light (v/c)
    speed_ratio = 0.99999987

    # Astronaut's lifespan and starting age
    lifespan = 150
    start_age = 22
    
    # The numerical options given in the question
    numerical_options = [72, 77, 81]

    # The distance to the Large Magellanic Cloud (LMC) is not specified.
    # We test a range of scientifically accepted values in light-years, as used by the candidate LLMs.
    distances_to_test = [159000, 160000, 163000]

    # --- Step 1: Calculate the Lorentz Factor (gamma) ---
    # This is a constant value based on the given speed.
    try:
        gamma = 1 / math.sqrt(1 - speed_ratio**2)
    except (ValueError, ZeroDivisionError) as e:
        return f"Incorrect: A fundamental calculation error occurred. Could not compute Lorentz factor. Error: {e}"

    # --- Step 2: Calculate Astronaut's Travel Time for Different Distances ---
    calculated_travel_times = []
    for distance_ly in distances_to_test:
        # Time in Earth's reference frame (t_earth)
        t_earth = distance_ly / speed_ratio
        
        # Time for the astronaut (proper time, t_astronaut)
        t_astronaut = t_earth / gamma
        calculated_travel_times.append(t_astronaut)

        # --- Step 3: Check Survival Constraint ---
        # The astronaut must survive the journey.
        remaining_lifespan = lifespan - start_age
        if t_astronaut >= remaining_lifespan:
            return (f"Incorrect: The 'Die' option would be correct for a distance of {distance_ly} ly. "
                    f"The journey takes {t_astronaut:.2f} years, which is not less than the remaining lifespan of {remaining_lifespan} years. "
                    f"This contradicts the final answer C.")

    # --- Step 4: Verify the Final Answer ---
    # The calculated times should consistently point to the given answer (81 years)
    # as the most plausible choice among the options.
    
    # Check if the provided answer is the closest option for all calculated times.
    for time in calculated_travel_times:
        # Find the closest option in the list to the calculated time
        closest_option = min(numerical_options, key=lambda x: abs(x - time))
        if closest_option != final_answer_value:
            return (f"Incorrect: For a distance of {distances_to_test[calculated_travel_times.index(time)]} ly, "
                    f"the calculated time is {time:.2f} years. The closest option is {closest_option} years, "
                    f"not the provided answer of {final_answer_value} years.")

    # If all checks pass, the answer is consistent and correct.
    return "Correct"

# Execute the check
result = check_correctness()
print(result)