import math

def check_answer():
    """
    Checks the correctness of the provided answer to the special relativity problem.
    """
    # Given values from the question
    v = 60000  # km/s
    t_proper = 50  # s
    c = 300000  # km/s, approximate speed of light

    # Options provided in the question
    options = {
        "A": 2880000,
        "B": 3060000,
        "C": 3000000,
        "D": 2940000
    }
    
    # The final answer provided by the LLM to be checked
    llm_answer_key = "B"

    # --- Step 1: Perform the physics calculation ---
    # This problem requires using the time dilation formula from special relativity.
    # The distance in the Galaxy's frame is d = v * t_galaxy.
    # The time in the Galaxy's frame (t_galaxy) is related to the proper time in the
    # ejecta's frame (t_proper) by the Lorentz factor (gamma).
    # t_galaxy = gamma * t_proper
    # gamma = 1 / sqrt(1 - (v/c)^2)

    # Calculate beta (v/c)
    beta = v / c
    
    # Check if velocity is less than the speed of light
    if beta >= 1:
        return "Calculation error: Velocity is greater than or equal to the speed of light."

    # Calculate the Lorentz factor (gamma)
    try:
        gamma = 1 / math.sqrt(1 - beta**2)
    except ValueError:
        return "Calculation error: Could not compute the square root, likely due to v >= c."

    # Calculate the time elapsed in the Galaxy's reference frame
    t_galaxy = gamma * t_proper

    # Calculate the distance traveled in the Galaxy's reference frame
    calculated_distance = v * t_galaxy

    # --- Step 2: Compare the calculated result with the options ---
    
    # Find the option closest to the calculated distance
    closest_option_key = min(options, key=lambda k: abs(options[k] - calculated_distance))
    closest_option_value = options[closest_option_key]

    # --- Step 3: Verify the LLM's answer ---
    
    # Check if the LLM's chosen answer matches the closest calculated option
    if llm_answer_key == closest_option_key:
        return "Correct"
    else:
        # Provide a detailed reason for the incorrectness
        reason = (
            f"The provided answer '{llm_answer_key}' is incorrect.\n"
            f"The calculation based on special relativity yields a distance of approximately {calculated_distance:,.2f} km.\n"
            f"The non-relativistic distance (ignoring time dilation) would be {v * t_proper:,.0f} km, which corresponds to option C.\n"
            f"Due to time dilation, the time in the Galaxy's frame is longer ({t_galaxy:.4f} s) than the proper time (50 s), so the distance must be greater than 3,000,000 km.\n"
            f"The calculated distance ({calculated_distance:,.2f} km) is closest to option '{closest_option_key}' ({closest_option_value:,.0f} km), not '{llm_answer_key}'."
        )
        return reason

# Run the check and print the result
print(check_answer())