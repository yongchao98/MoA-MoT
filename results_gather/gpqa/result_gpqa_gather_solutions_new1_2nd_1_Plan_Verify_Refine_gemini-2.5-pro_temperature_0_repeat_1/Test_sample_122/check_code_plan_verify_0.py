import math

def check_supernova_distance():
    """
    This function checks the correctness of the answer to the special relativity problem.

    It calculates the theoretical distance traveled by the ejecta in the Galaxy's
    reference frame and compares it to the provided options to verify the chosen answer.
    """
    # --- Problem Parameters ---
    # Velocity of the ejecta (v) in km/s
    v = 60000
    # Proper time in the ejecta's frame (t0) in seconds
    t0 = 50
    # Speed of light (c) in km/s
    c = 300000

    # --- Options from the Question ---
    # The options are provided as a dictionary for easy lookup.
    options = {
        "A": 2940000,
        "B": 3060000,
        "C": 3000000,
        "D": 2880000
    }

    # The final answer provided by the LLM to be checked.
    llm_answer_key = "B"

    # --- Physics Calculation ---
    # Step 1: Calculate the non-relativistic (classical) distance for a sanity check.
    classical_distance = v * t0
    if options["C"] != classical_distance:
        return f"Constraint check failed: Option C should represent the classical distance ({classical_distance:,} km), but it is {options['C']:,} km."

    # Step 2: Calculate the Lorentz factor (gamma).
    # beta is the ratio of the object's velocity to the speed of light.
    beta = v / c
    # gamma is the factor by which time, length, and relativistic mass change for a moving object.
    gamma = 1 / math.sqrt(1 - beta**2)

    # Step 3: Calculate the time elapsed in the Galaxy's frame (dilated time).
    t_galaxy = gamma * t0

    # Step 4: Calculate the distance traveled in the Galaxy's frame.
    d_galaxy_calculated = v * t_galaxy

    # --- Verification ---
    # Step 5: Check the logical sanity constraint.
    # Due to time dilation, the time in the galaxy frame is longer, so the distance
    # must be greater than the classical distance.
    if options[llm_answer_key] <= classical_distance:
        return (f"The chosen answer '{llm_answer_key}' ({options[llm_answer_key]:,} km) is incorrect. "
                f"It violates the principle of time dilation, which requires the distance to be "
                f"greater than the non-relativistic distance of {classical_distance:,} km.")

    # Step 6: Find which option is numerically closest to the calculated result.
    # We calculate the absolute difference between our calculated distance and each option's value.
    # The key with the minimum difference is the correct option.
    closest_option_key = min(options, key=lambda k: abs(options[k] - d_galaxy_calculated))

    # Step 7: Compare the LLM's answer with the closest calculated option.
    if llm_answer_key == closest_option_key:
        return "Correct"
    else:
        return (f"The answer '{llm_answer_key}' is incorrect. "
                f"The calculated relativistic distance is approximately {d_galaxy_calculated:,.2f} km. "
                f"The closest option is '{closest_option_key}' ({options[closest_option_key]:,} km), "
                f"not '{llm_answer_key}' ({options[llm_answer_key]:,} km).")

# Run the check and print the result.
result = check_supernova_distance()
print(result)