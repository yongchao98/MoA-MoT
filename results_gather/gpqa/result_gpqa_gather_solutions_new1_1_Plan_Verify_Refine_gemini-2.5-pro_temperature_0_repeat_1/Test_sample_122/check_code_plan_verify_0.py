import math

def check_correctness():
    """
    Checks the correctness of the final answer for the special relativity problem.

    The problem involves calculating the distance traveled by a supernova ejecta
    in the Galaxy's reference frame, given the time elapsed in the ejecta's frame.
    This requires using the time dilation formula from special relativity.

    Steps:
    1. Define the given physical quantities.
    2. Define the multiple-choice options.
    3. Calculate the Lorentz factor (gamma).
    4. Calculate the time elapsed in the Galaxy's frame (dilated time).
    5. Calculate the distance traveled in the Galaxy's frame.
    6. Compare the calculated distance to the options to find the closest one.
    7. Check if the closest option matches the provided answer.
    """
    # --- Problem Constraints and Given Values ---
    v = 60000  # km/s (relative velocity)
    t0 = 50      # s (proper time in the ejecta's frame)
    c = 300000   # km/s (approximate speed of light)

    # --- Options from the Question ---
    # The question provides these options:
    # A) 3 000 000 km
    # B) 3 060 000 km
    # C) 2 880 000 km
    # D) 2 940 000 km
    options = {
        "A": 3000000,
        "B": 3060000,
        "C": 2880000,
        "D": 2940000
    }

    # The final answer provided by the LLM to be checked
    llm_answer_letter = "B"

    # --- Physics Calculation ---
    # The distance in the Galaxy's frame is d = v * t_galaxy.
    # We need to find t_galaxy, the time elapsed in the Galaxy's frame.
    # This is found using the time dilation formula: t_galaxy = gamma * t0.

    # 1. Calculate the Lorentz factor (gamma)
    try:
        beta_squared = (v / c)**2
        if beta_squared >= 1:
            return "Calculation error: Velocity cannot be equal to or greater than the speed of light."
        gamma = 1 / math.sqrt(1 - beta_squared)
    except Exception as e:
        return f"An error occurred during calculation: {e}"

    # 2. Calculate the time in the Galaxy's frame (dilated time)
    t_galaxy = gamma * t0

    # 3. Calculate the distance in the Galaxy's frame
    calculated_distance = v * t_galaxy

    # --- Verification ---
    # Find which option is numerically closest to our calculated distance.
    closest_option = None
    min_difference = float('inf')

    for letter, value in options.items():
        difference = abs(calculated_distance - value)
        if difference < min_difference:
            min_difference = difference
            closest_option = letter

    # Check if the LLM's answer matches the closest option.
    # We also check if the difference is reasonably small (e.g., less than 1% of the value),
    # which accounts for potential rounding in the problem's options.
    if closest_option == llm_answer_letter:
        tolerance = 0.01 * options[llm_answer_letter] # 1% tolerance
        if min_difference <= tolerance:
            return "Correct"
        else:
            return (f"The provided answer '{llm_answer_letter}' is the closest option, but the numerical difference is large. "
                    f"Calculated distance: {calculated_distance:.2f} km. "
                    f"Option {llm_answer_letter}: {options[llm_answer_letter]} km. "
                    f"Difference: {min_difference:.2f} km.")
    else:
        return (f"Incorrect. The provided answer is '{llm_answer_letter}', but the calculation shows '{closest_option}' is the correct choice. "
                f"The calculated distance is {calculated_distance:.2f} km, which is closest to option {closest_option} ({options[closest_option]} km). "
                f"The non-relativistic distance (v*t0) would be {v*t0} km (Option A), which is incorrect.")

# Execute the check and print the result
result = check_correctness()
print(result)