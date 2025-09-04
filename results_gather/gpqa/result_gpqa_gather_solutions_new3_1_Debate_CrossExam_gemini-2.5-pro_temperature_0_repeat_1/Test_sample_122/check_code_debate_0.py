import math

def check_correctness():
    """
    This function checks the correctness of the provided answer to the special relativity problem.

    The problem involves calculating the distance traveled by a supernova ejecta in the Galaxy's
    reference frame, given the time elapsed in the ejecta's own frame. This requires using
    the time dilation formula from special relativity.

    Steps:
    1. Define the given physical constants and values.
    2. Calculate the Lorentz factor (gamma).
    3. Calculate the dilated time as measured in the Galaxy's frame.
    4. Calculate the distance traveled in the Galaxy's frame.
    5. Compare the calculated distance to the provided options to find the closest match.
    6. Check if the provided answer corresponds to this closest match.
    """
    # Given values from the question
    v = 60000  # Relative velocity in km/s
    t_proper = 50  # Proper time in the ejecta's frame in seconds
    c = 300000  # Speed of light in km/s (standard approximation for such problems)

    # The options provided in the question
    options = {
        "A": 2880000,
        "B": 3000000,
        "C": 2940000,
        "D": 3060000
    }

    # The final answer from the LLM to be checked
    llm_answer_choice = "D"

    # --- Physics Calculation ---

    # 1. Calculate the Lorentz factor (gamma)
    # gamma = 1 / sqrt(1 - v^2/c^2)
    try:
        beta_squared = (v / c) ** 2
        if beta_squared >= 1:
            return "Calculation error: Velocity is greater than or equal to the speed of light."
        gamma = 1 / math.sqrt(1 - beta_squared)
    except Exception as e:
        return f"An error occurred during calculation: {e}"

    # 2. Calculate the time elapsed in the Galaxy's frame (dilated time)
    # t_galaxy = gamma * t_proper
    t_galaxy = gamma * t_proper

    # 3. Calculate the distance traveled in the Galaxy's frame
    # distance = velocity * time_in_galaxy_frame
    calculated_distance = v * t_galaxy

    # --- Verification ---

    # Find the option that is numerically closest to our calculated distance
    closest_option_key = min(options, key=lambda k: abs(options[k] - calculated_distance))

    # Check if the LLM's answer matches the correct calculation
    if llm_answer_choice == closest_option_key:
        # As a final check, ensure the chosen answer is not the non-relativistic one by mistake.
        non_relativistic_distance = v * t_proper
        if abs(options[llm_answer_choice] - non_relativistic_distance) < 1: # Using a small tolerance for float comparison
            return (f"Incorrect. The chosen answer {llm_answer_choice} ({options[llm_answer_choice]} km) corresponds to the "
                    f"non-relativistic calculation (v * t_proper = {non_relativistic_distance} km). "
                    f"The problem requires applying time dilation, which leads to a distance of approximately {calculated_distance:.0f} km.")
        return "Correct"
    else:
        return (f"Incorrect. The calculation shows the distance is approximately {calculated_distance:.0f} km. "
                f"The closest option is {closest_option_key} ({options[closest_option_key]} km), "
                f"but the provided answer was {llm_answer_choice} ({options[llm_answer_choice]} km).")

# Run the check and print the result
result = check_correctness()
print(result)