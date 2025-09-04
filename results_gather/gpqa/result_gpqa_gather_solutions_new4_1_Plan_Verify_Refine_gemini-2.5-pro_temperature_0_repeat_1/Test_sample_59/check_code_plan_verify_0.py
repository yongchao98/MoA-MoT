import math

def check_correctness():
    """
    This function checks the correctness of the provided answer to the special relativity problem.

    The core of the problem is calculating the time experienced by the astronaut (proper time)
    using the time dilation formula. A key piece of information, the distance to the
    Large Magellanic Cloud (LMC), is not provided in the question and must be estimated
    from standard astronomical data.
    """

    # --- Problem Constants & Given Data ---
    v_ratio = 0.99999987  # Speed of the spacecraft as a fraction of the speed of light (v/c)
    age_start = 22        # Astronaut's starting age in years
    lifespan = 150        # Alien's average lifespan in years

    # --- LLM's Answer to be Checked ---
    # The provided answer is 'B', which corresponds to a travel time of 81 years.
    llm_answer_time = 81

    # --- Verification Steps ---

    # 1. Handle the missing distance information.
    # The distance to the LMC is known to be approximately 158,000 to 163,000 light-years.
    # A common approximation is 160,000 light-years. We will use this value to verify the calculation,
    # as was done in the provided reasoning.
    distance_ly = 160000

    # 2. Perform the physics calculation for time dilation.
    # Formula for proper time (t_astronaut): t_astronaut = (d / v) * sqrt(1 - (v/c)^2)
    # This can be calculated in two steps:
    # a) Calculate time in Earth's frame (t_earth).
    # b) Apply the Lorentz factor to find the astronaut's time.
    try:
        # Time in Earth's frame (in years)
        t_earth = distance_ly / v_ratio
        
        # The time dilation factor
        dilation_factor = math.sqrt(1 - v_ratio**2)
        
        # Time for the astronaut (proper time)
        calculated_t_astronaut = t_earth * dilation_factor
    except (ValueError, ZeroDivisionError) as e:
        return f"An error occurred during the physics calculation: {e}"

    # 3. Check if the LLM's answer is the most plausible option.
    # We compare our calculated time to the available numerical options to see which is closest.
    numerical_options = [72, 77, 81]
    closest_option = min(numerical_options, key=lambda x: abs(x - calculated_t_astronaut))

    if closest_option != llm_answer_time:
        return (f"Incorrect. Using a standard distance of {distance_ly} light-years, the calculated travel time "
                f"is {calculated_t_astronaut:.2f} years. The closest option to this value is {closest_option} years, "
                f"not the provided answer of {llm_answer_time} years.")

    # 4. Check the survival constraint.
    # The astronaut must survive the trip for the answer to be valid.
    age_on_arrival = age_start + llm_answer_time
    if age_on_arrival >= lifespan:
        return (f"Incorrect. The answer of {llm_answer_time} years implies the astronaut would arrive at age "
                f"{age_on_arrival}, which is not within their {lifespan}-year lifespan. This contradicts the premise "
                f"that the journey is completed and a time is given.")

    # If all checks pass, the answer is consistent and correct.
    return "Correct"

# The final output of the check.
print(check_correctness())