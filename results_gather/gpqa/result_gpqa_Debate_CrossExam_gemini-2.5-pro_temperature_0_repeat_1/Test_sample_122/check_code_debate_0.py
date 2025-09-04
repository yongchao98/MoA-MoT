import math

def check_answer():
    """
    This function checks the correctness of the provided answer to the special relativity problem.
    """
    # --- Define constants and given values from the problem ---
    # Velocity of the ejecta relative to the Galaxy (v)
    v = 60000  # in km/s
    # Time passed in the ejecta's reference frame (proper time, Δt')
    dt_prime = 50  # in seconds
    # Speed of light (c)
    c = 300000  # in km/s

    # --- Define the options and the provided answer ---
    options = {
        "A": 2880000,
        "B": 3000000,
        "C": 3060000,
        "D": 2940000
    }
    llm_answer_key = "C"

    # --- Physics Calculation ---
    # The question asks for the distance traveled by the ejecta in the Galaxy's reference frame.
    # Let S be the Galaxy's frame (stationary) and S' be the ejecta's frame (moving).
    # The distance 'd' in frame S is given by d = v * Δt, where Δt is the time elapsed in frame S.
    # The time Δt in frame S is related to the proper time Δt' in frame S' by the time dilation formula:
    # Δt = γ * Δt'
    # where γ (the Lorentz factor) is 1 / sqrt(1 - v²/c²)

    # Combining these, the formula for the distance is:
    # d = v * (γ * Δt') = v * Δt' / sqrt(1 - v²/c²)

    try:
        # Calculate the term inside the square root
        beta_squared = (v / c)**2
        if beta_squared >= 1:
            return "Calculation error: The velocity is greater than or equal to the speed of light."
        
        denominator = math.sqrt(1 - beta_squared)
        
        # Calculate the distance
        calculated_distance = (v * dt_prime) / denominator

    except Exception as e:
        return f"An error occurred during calculation: {e}"

    # --- Verification ---
    # Find which of the given options is closest to our calculated distance.
    closest_option_key = None
    min_difference = float('inf')

    for key, value in options.items():
        difference = abs(calculated_distance - value)
        if difference < min_difference:
            min_difference = difference
            closest_option_key = key

    # Check if the LLM's answer matches the closest option.
    if closest_option_key == llm_answer_key:
        return "Correct"
    else:
        llm_answer_value = options[llm_answer_key]
        closest_option_value = options[closest_option_key]
        return (f"Incorrect. The provided answer is {llm_answer_key} ({llm_answer_value:,} km), but the calculation shows the closest answer is {closest_option_key} ({closest_option_value:,} km).\n"
                f"Detailed calculation: The distance traveled in the Galaxy's frame is d = v * Δt' / sqrt(1 - v²/c²).\n"
                f"d = (60,000 km/s * 50 s) / sqrt(1 - (60,000/300,000)²)\n"
                f"d = 3,000,000 / sqrt(1 - 0.2²)\n"
                f"d = 3,000,000 / sqrt(0.96)\n"
                f"d ≈ {calculated_distance:,.2f} km.\n"
                f"This calculated value is closest to option {closest_option_key}.")

# Run the check and print the result
result = check_answer()
print(result)