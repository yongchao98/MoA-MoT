import math

def check_correctness():
    """
    Checks the correctness of the final answer to the physics problem.

    The final answer concludes that the journey takes approximately 81 years
    from the astronaut's perspective, and the astronaut survives. This check
    verifies these two key points.
    """

    # --- Given constants and parameters ---
    v_ratio = 0.99999987  # v/c
    initial_age = 22
    lifespan = 150

    # --- Options from the question (numerical values) ---
    # The options are 72, 77, 81 years, or death.
    # The answer to check is that 81 years is the correct choice.
    options = {'72 years': 72, '77 years': 77, '81 years': 81}
    target_answer_value = 81

    # --- Handle ambiguity in distance ---
    # The distance to the Large Magellanic Cloud (LMC) is not given.
    # We test a range of commonly accepted astronomical values.
    distances_ly = [159000, 160000, 163000]

    # --- Perform the physics calculation ---
    try:
        # The Lorentz factor is constant as it only depends on velocity.
        lorentz_factor = 1 / math.sqrt(1 - v_ratio**2)
    except ValueError:
        return "Calculation Error: The velocity ratio v/c cannot be 1 or greater."

    for distance in distances_ly:
        # Calculate time in Earth's reference frame
        delta_t_earth = distance / v_ratio

        # Calculate time in the astronaut's reference frame (proper time)
        delta_t_astronaut = delta_t_earth / lorentz_factor

        # --- Check Constraint 1: Survival ---
        final_age = initial_age + delta_t_astronaut
        if final_age >= lifespan:
            return (f"Incorrect. The answer states the astronaut survives, but this is not guaranteed. "
                    f"With a distance of {distance} ly, the journey takes {delta_t_astronaut:.2f} years, "
                    f"making the astronaut's final age {final_age:.2f}, which exceeds the {lifespan}-year lifespan.")

        # --- Check Constraint 2: Closest Option ---
        # Find which of the numerical options is closest to our calculated time.
        closest_option_value = min(options.values(), key=lambda val: abs(val - delta_t_astronaut))

        if closest_option_value != target_answer_value:
            return (f"Incorrect. The answer claims 81 years is the correct choice, but this is not always the case. "
                    f"With a distance of {distance} ly, the calculated time is {delta_t_astronaut:.2f} years. "
                    f"The closest option to this value is {closest_option_value} years, not {target_answer_value} years.")

    # If the code reaches this point, it means for all tested plausible distances:
    # 1. The astronaut survives the journey.
    # 2. The calculated travel time is always closest to 81 years.
    # Therefore, the logic and conclusion of the provided answer are sound.
    return "Correct"

# Run the check
result = check_correctness()
print(result)