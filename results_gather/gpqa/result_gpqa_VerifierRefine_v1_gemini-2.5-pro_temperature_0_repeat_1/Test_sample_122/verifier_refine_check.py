import math

def check_relativity_answer():
    """
    This function checks the correctness of the answer to the special relativity problem.

    It calculates the distance traveled by the supernova ejecta in the Galaxy's
    reference frame using the principles of time dilation.
    """
    # --- Given parameters from the question ---
    # Relative velocity (v) in km/s
    v = 60000.0
    # Proper time (Δt₀) in seconds, measured in the ejecta's frame
    proper_time_t0 = 50.0
    # Speed of light (c) in km/s
    c = 300000.0

    # --- The answer provided by the LLM ---
    # The LLM chose option C
    llm_answer_choice = 'C'
    options = {
        'A': 2940000.0,
        'B': 2880000.0,
        'C': 3060000.0,
        'D': 3000000.0
    }
    llm_answer_value = options.get(llm_answer_choice)

    if llm_answer_value is None:
        return f"Invalid answer choice '{llm_answer_choice}' provided."

    # --- Physics Calculation ---

    # 1. Calculate the Lorentz factor (gamma)
    # The formula is γ = 1 / sqrt(1 - v²/c²)
    try:
        beta_squared = (v / c)**2
        gamma = 1 / math.sqrt(1 - beta_squared)
    except ValueError:
        return "Calculation Error: The velocity 'v' cannot be greater than or equal to the speed of light 'c'."

    # 2. Calculate the time elapsed in the Galaxy's frame (Δt) using time dilation
    # The formula is Δt = γ * Δt₀
    galaxy_time_t = gamma * proper_time_t0

    # 3. Calculate the distance traveled in the Galaxy's frame (d)
    # The formula is d = v * Δt
    calculated_distance = v * galaxy_time_t

    # --- Verification ---
    # Find which option is numerically closest to our calculated result.
    closest_option = min(options, key=lambda k: abs(options[k] - calculated_distance))

    # Check if the LLM's chosen answer is the closest option.
    if closest_option == llm_answer_choice:
        # The answer is correct as it represents the closest numerical approximation.
        # The exact calculated value is ~3,061,862 km, which rounds well to 3,060,000 km.
        return "Correct"
    else:
        return (f"Incorrect. The calculation shows the distance should be approximately {calculated_distance:,.2f} km. "
                f"The closest option is '{closest_option}' ({options[closest_option]:,.0f} km), "
                f"but the provided answer was '{llm_answer_choice}' ({llm_answer_value:,.0f} km).")

# Execute the check and print the result
result = check_relativity_answer()
print(result)