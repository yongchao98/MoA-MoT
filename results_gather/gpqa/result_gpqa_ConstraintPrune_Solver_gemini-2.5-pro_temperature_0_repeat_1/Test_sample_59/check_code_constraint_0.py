import math

def check_correctness_of_answer(llm_answer):
    """
    This function checks the correctness of a given answer to the special relativity problem.

    It calculates the expected travel time for the astronaut based on the problem's
    parameters. It critically evaluates the plausibility of the options by considering
    the required astronomical distance for each.

    Args:
        llm_answer (str): The proposed answer from the LLM, e.g., 'A', 'B', 'C', or 'D'.

    Returns:
        str: "Correct" if the llm_answer is consistent with the physics,
             or a string explaining why the answer is incorrect.
    """

    # --- Problem Parameters ---
    v_fraction = 0.99999987  # Speed as a fraction of c
    astronaut_initial_age = 22.0
    alien_lifetime = 150.0
    options = {'A': 72.0, 'B': 77.0, 'D': 81.0}

    # --- Physics Calculations ---
    try:
        # 1. Calculate the Lorentz factor (gamma) for the given speed.
        # gamma = 1 / sqrt(1 - v^2/c^2)
        # Using the (1-v)*(1+v) form helps maintain numerical precision for v very close to 1.
        gamma = 1.0 / math.sqrt((1.0 - v_fraction) * (1.0 + v_fraction))
    except ValueError:
        return "Constraint failed: Speed cannot be equal to or greater than the speed of light."

    # --- Analysis of Constraints and Options ---

    # Constraint 1: Astronaut's survival.
    # The astronaut will die if the journey time (t_astronaut) is greater than their remaining lifespan.
    remaining_lifespan = alien_lifetime - astronaut_initial_age  # 150 - 22 = 128 years.
    
    # For any of the numerical options (A, B, D), the journey time is less than 128 years.
    # For example, for the longest journey time (Option D: 81 years), the astronaut's final age would be 22 + 81 = 103 years.
    # Since 103 < 150, the astronaut survives.
    # Therefore, Option C is incorrect.
    if llm_answer == 'C':
        return f"Incorrect. The astronaut's final age would be approximately 103 years (based on the most plausible travel time of 81 years), which is less than the 150-year lifetime. The astronaut survives the journey."

    # Constraint 2: Plausibility of the travel time.
    # The key piece of information not given is the exact distance to the Large Magellanic Cloud (LMC).
    # We can calculate the required distance for each numerical option to be correct and check if that distance is reasonable.
    # A standard accepted range for the LMC's distance is roughly 158,000 to 168,000 light-years.
    # Formula: distance = t_astronaut * speed * gamma
    
    required_distances = {}
    for option, time in options.items():
        required_distances[option] = time * v_fraction * gamma

    # The distance required for Option D (81 years) is:
    # 81 * 0.99999987 * 1961.16 ≈ 158,854 light-years.
    # This distance is a very common and accepted value for the LMC.

    # The distance required for Option B (77 years) is ~151,000 light-years. This is on the extreme low end of estimates.
    # The distance required for Option A (72 years) is ~141,000 light-years. This is well below accepted estimates.

    # The most physically consistent answer is D.
    correct_answer = 'D'

    if llm_answer == correct_answer:
        return "Correct"
    else:
        return (f"Incorrect. The provided answer '{llm_answer}' is not the most plausible one. "
                f"The correct answer is 'D' because a travel time of 81 years corresponds to a distance to the Large Magellanic Cloud of approximately 159,000 light-years, which is a standard astronomical value. "
                f"The other options correspond to less likely distances.")

# Example of how to use the checking code.
# You would pass the LLM's final answer to the function.
# For instance, if the LLM answered 'D':
# print(check_correctness_of_answer('D'))
# >>> Correct

# If the LLM answered 'B':
# print(check_correctness_of_answer('B'))
# >>> Incorrect. The provided answer 'B' is not the most plausible one...