import math

def check_correctness():
    """
    Checks the correctness of the LLM's answer to the special relativity problem.
    """
    # 1. Define constants and given values from the question
    v_ratio = 0.99999987  # v/c, speed of the spacecraft as a fraction of the speed of light
    astronaut_initial_age = 22  # years
    alien_lifespan = 150  # years

    # The distance to the Large Magellanic Cloud (LMC) is an unstated but necessary value.
    # Standard values range from 159,000 to 163,000 light-years.
    # The provided solution correctly infers that ~160,000 light-years is the intended value.
    distance_ly = 160000.0  # in light-years

    # The multiple-choice options from the question
    options = {'A': 72, 'B': 'The astronaut will die before reaching to the Earth.', 'C': 77, 'D': 81}
    
    # The final answer provided by the LLM to be checked
    llm_answer_key = 'D'

    # 2. Perform the physics calculation
    # Step 2a: Calculate the Lorentz factor (gamma)
    # gamma = 1 / sqrt(1 - (v/c)^2)
    try:
        lorentz_factor = 1 / math.sqrt(1 - v_ratio**2)
    except ValueError:
        return "Calculation Error: Velocity ratio implies speed >= c, which is physically impossible."

    # Step 2b: Calculate the travel time in Earth's reference frame (delta_t)
    # Since distance is in light-years and speed is a fraction of c, the time in years is d/v_ratio.
    time_earth_frame = distance_ly / v_ratio

    # Step 2c: Calculate the travel time in the astronaut's frame (proper time, delta_t0)
    time_astronaut_frame = time_earth_frame / lorentz_factor

    # 3. Check the constraints and correctness of the answer
    # Step 3a: Check the survival constraint
    astronaut_final_age = astronaut_initial_age + time_astronaut_frame
    astronaut_survives = astronaut_final_age < alien_lifespan

    # Step 3b: Evaluate the LLM's answer based on our calculation
    # The LLM chose 'D', which implies the astronaut survives. Let's verify this.
    if not astronaut_survives:
        return (f"Incorrect. The LLM's answer implies the astronaut survives. However, the calculation shows "
                f"the journey takes {time_astronaut_frame:.2f} years, making the astronaut's final age "
                f"{astronaut_final_age:.2f} years, which exceeds their 150-year lifespan. The correct answer "
                f"should have been B.")

    # Step 3c: Find the closest numerical option to our calculated time
    closest_option_key = None
    min_difference = float('inf')

    for key, value in options.items():
        if isinstance(value, (int, float)):
            difference = abs(time_astronaut_frame - value)
            if difference < min_difference:
                min_difference = difference
                closest_option_key = key
    
    # Step 3d: Final verification
    if llm_answer_key == closest_option_key:
        # The LLM chose the option that is mathematically closest to the calculated result.
        return "Correct"
    else:
        return (f"Incorrect. The calculated travel time for the astronaut is approximately {time_astronaut_frame:.2f} years. "
                f"Based on the options, the closest answer is {options[closest_option_key]} years (Option {closest_option_key}). "
                f"The provided answer was Option {llm_answer_key} ({options[llm_answer_key]} years).")

# Run the check and print the result
result = check_correctness()
print(result)