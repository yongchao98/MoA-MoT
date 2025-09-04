import math

def check_relativity_problem():
    """
    This function verifies the solution to the special relativity problem.

    It calculates:
    1. The Lorentz factor based on the given velocity.
    2. The travel time from Earth's perspective, assuming a standard distance to the LMC.
    3. The travel time from the astronaut's perspective (proper time).
    4. The astronaut's final age to check the survival constraint.
    5. Compares the calculated result with the provided answer.
    """
    # --- Problem Parameters & Constraints ---
    v_over_c = 0.99999987  # Speed as a fraction of the speed of light
    astronaut_initial_age = 22  # years
    alien_lifespan = 150  # years

    # The distance to the Large Magellanic Cloud (LMC) is not given.
    # We use a standard, commonly cited value of 160,000 light-years,
    # as this value aligns best with the provided options.
    distance_ly = 160000.0

    # The final answer provided is 'A', which corresponds to 81 years in the question's option list.
    options = {
        "A": 81.0,
        "B": 72.0,
        "C": 77.0,
        "D": "The astronaut will die before reaching to the Earth."
    }
    provided_answer_key = "A"
    expected_value = options[provided_answer_key]

    # --- Step 1: Calculate the Lorentz Factor (gamma) ---
    # gamma = 1 / sqrt(1 - v^2/c^2)
    try:
        gamma = 1 / math.sqrt(1 - v_over_c**2)
    except ValueError:
        return "Calculation Error: Velocity cannot be equal to or greater than the speed of light."

    # --- Step 2: Calculate Travel Time in Earth's Reference Frame (delta_t) ---
    # delta_t = distance / velocity
    # Since distance is in light-years and velocity is a fraction of c, the time is in years.
    delta_t_earth = distance_ly / v_over_c

    # --- Step 3: Calculate Travel Time for the Astronaut (Proper Time, delta_t_0) ---
    # delta_t_0 = delta_t / gamma
    delta_t_astronaut = delta_t_earth / gamma

    # --- Step 4: Check the Survival Constraint ---
    final_age = astronaut_initial_age + delta_t_astronaut
    survives = final_age < alien_lifespan

    # --- Step 5: Verify the Answer ---
    # First, check if the survival prediction is correct.
    # The provided answer 'A' (81 years) implies the astronaut survives.
    if not survives:
        if expected_value == options["D"]:
            return "Correct"
        else:
            return (f"Incorrect. The calculation shows the astronaut's final age would be {final_age:.2f} years, "
                    f"which exceeds their lifespan of {alien_lifespan} years. The astronaut would not survive. "
                    f"The provided answer implies survival.")

    # The astronaut survives, so the answer should be a numerical value.
    if expected_value == options["D"]:
        return (f"Incorrect. The calculation shows the astronaut's final age would be {final_age:.2f} years, "
                f"which is within their lifespan of {alien_lifespan} years. The astronaut would survive. "
                f"The provided answer incorrectly states the astronaut would die.")

    # Find the closest numerical option to our calculated time.
    numerical_options = {k: v for k, v in options.items() if isinstance(v, (int, float))}
    closest_option_key = min(numerical_options, key=lambda k: abs(numerical_options[k] - delta_t_astronaut))

    # Check if the provided answer is the closest option.
    if provided_answer_key != closest_option_key:
        return (f"Incorrect. The calculated travel time for the astronaut is {delta_t_astronaut:.2f} years. "
                f"This is closest to option {closest_option_key} ({options[closest_option_key]} years), "
                f"not the provided answer {provided_answer_key} ({expected_value} years).")

    # The provided answer is the closest. Check if the difference is reasonably small.
    # A small tolerance accounts for rounding or using a slightly different distance to LMC.
    if abs(delta_t_astronaut - expected_value) > 2.0:
        return (f"Incorrect. While the provided answer '{expected_value} years' is the closest option, "
                f"the calculated time is {delta_t_astronaut:.2f} years. The difference is larger than expected, "
                f"which might indicate an issue with the problem's premises or options.")

    # All checks passed.
    return "Correct"

# Run the check and print the result.
result = check_relativity_problem()
print(result)