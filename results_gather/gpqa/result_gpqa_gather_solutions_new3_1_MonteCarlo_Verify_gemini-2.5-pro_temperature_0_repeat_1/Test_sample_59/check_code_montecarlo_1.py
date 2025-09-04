import math

def check_correctness():
    """
    Checks the correctness of the provided answer to the physics problem.
    The problem asks for the travel time from the astronaut's perspective.
    The provided answer to check is 'A', which corresponds to 81 years.
    """

    # --- Define problem parameters and constraints ---
    v_c_ratio = 0.99999987  # Speed as a fraction of the speed of light
    astronaut_age_start = 22
    alien_lifespan = 150
    
    # The options given in the question are:
    # A) 81 years
    # B) The astronaut will die before reaching to the Earth.
    # C) 72 years
    # D) 77 years
    
    # The answer to be checked is 'A', which corresponds to a travel time of 81 years.
    answer_to_check_value = 81.0

    # --- Handle ambiguity in distance ---
    # The distance to the Large Magellanic Cloud (LMC) is not specified in the question.
    # Astronomical sources give values between ~159,000 and ~163,000 light-years.
    # A common approach for such problems is to assume the intended distance is one that
    # leads cleanly to one of the options. We will test with a common value of 160,000 light-years,
    # which is used in several of the candidate LLM answers.
    distance_in_ly = 160000.0

    # --- Perform the physics calculation ---
    # 1. Calculate the Lorentz factor (gamma)
    try:
        gamma = 1.0 / math.sqrt(1 - v_c_ratio**2)
    except ValueError:
        return "Calculation Error: The speed ratio results in a math domain error."

    # 2. Calculate time in Earth's reference frame (delta_t)
    delta_t_earth = distance_in_ly / v_c_ratio

    # 3. Calculate time in the astronaut's reference frame (proper time, delta_t0)
    delta_t_astronaut = delta_t_earth / gamma

    # --- Verify the answer against constraints and calculations ---
    # Constraint 1: Survival
    # The astronaut's age upon arrival would be:
    age_at_arrival = astronaut_age_start + delta_t_astronaut
    
    if age_at_arrival >= alien_lifespan:
        # This would mean option B is correct.
        # The answer to check is 'A' (81 years), which implies survival.
        # If our calculation shows the astronaut dies, then 'A' is incorrect.
        return (f"Incorrect. The calculated travel time is {delta_t_astronaut:.2f} years. "
                f"The astronaut's age on arrival would be {age_at_arrival:.2f}, which "
                f"exceeds the lifespan of {alien_lifespan}. Therefore, option B ('The astronaut will die') "
                f"would be correct, not A.")

    # Constraint 2: Closeness to options
    # The astronaut survives, so option B is incorrect. Now we check which numerical option is closest.
    # The numerical options are 81, 72, 77.
    numerical_options = [81.0, 72.0, 77.0]
    
    # Find the option value that is closest to our calculated time
    closest_option = min(numerical_options, key=lambda x: abs(x - delta_t_astronaut))

    # The answer to check is 81 years. Our calculation gives ~81.58 years.
    # The closest numerical option to 81.58 is 81.
    if closest_option == answer_to_check_value:
        # The calculation confirms that 81 years is the most plausible answer among the choices.
        # Since the provided answer is 'A' (81 years), it is correct.
        return "Correct"
    else:
        # This case would trigger if, for example, the calculation resulted in 76 years,
        # making 77 the closest option.
        return (f"Incorrect. Using a distance of {distance_in_ly} ly, the calculated travel time is "
                f"{delta_t_astronaut:.2f} years. This is closest to the option '{closest_option} years', "
                f"not the provided answer of '{answer_to_check_value} years'.")

# Execute the check and print the result.
result = check_correctness()
print(result)