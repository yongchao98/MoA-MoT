import math

def check_answer():
    """
    Checks the correctness of the answer to the special relativity problem.
    """
    # --- Define constants and given values from the question ---
    # Velocity of the ejecta relative to the Galaxy
    v = 60000  # km/s
    # Proper time elapsed in the ejecta's reference frame
    dt0 = 50   # s
    # Speed of light
    c = 300000 # km/s

    # --- The final answer provided by the LLM to be checked ---
    # The LLM's final answer is 'B', which corresponds to 3,060,000 km.
    llm_answer_letter = 'B'
    
    # --- The options as presented in the question ---
    options = {
        'A': 2880000,
        'B': 3060000,
        'C': 2940000,
        'D': 3000000
    }

    # --- Step 1: Perform the physics calculation ---
    # The distance in the Galaxy's frame is d = v * Δt, where Δt is the time in the Galaxy's frame.
    # Δt is related to the proper time Δt₀ by the time dilation formula: Δt = γ * Δt₀
    # where γ (gamma) is the Lorentz factor: γ = 1 / sqrt(1 - v²/c²)
    
    # Calculate the Lorentz factor (gamma)
    try:
        beta_squared = (v / c)**2
        if beta_squared >= 1:
            return "Calculation error: Velocity cannot be equal to or greater than the speed of light."
        gamma = 1 / math.sqrt(1 - beta_squared)
    except Exception as e:
        return f"An error occurred during calculation: {e}"

    # Calculate the distance traveled in the Galaxy's reference frame
    calculated_distance = v * gamma * dt0

    # --- Step 2: Verify the correctness of the provided answer ---
    
    # Check if the provided answer letter is valid
    if llm_answer_letter not in options:
        return f"Invalid answer format. The answer '{llm_answer_letter}' is not one of the options A, B, C, D."

    llm_answer_value = options[llm_answer_letter]

    # Constraint 1: Sanity check based on physics principles.
    # The classical (non-relativistic) distance would be v * dt0.
    # Due to time dilation, the time in the Galaxy's frame is longer, so the distance must be greater than the classical distance.
    classical_distance = v * dt0
    if llm_answer_value <= classical_distance:
        return (f"Incorrect. The chosen answer {llm_answer_value} km is not greater than the classical distance of {classical_distance} km. "
                f"Special relativity dictates the distance must be greater due to time dilation.")

    # Constraint 2: The chosen option should be the numerically closest to the precise calculation.
    # Find the option that is closest to our calculated distance.
    closest_option_letter = min(options, key=lambda k: abs(options[k] - calculated_distance))

    # Check if the LLM's answer matches the closest option.
    if llm_answer_letter == closest_option_letter:
        # The answer is correct as it's the closest option and satisfies the physical constraints.
        return "Correct"
    else:
        # The answer is incorrect because it's not the closest option.
        return (f"Incorrect. The precise calculated distance is approximately {calculated_distance:.2f} km. "
                f"The closest option to this value is {closest_option_letter} ({options[closest_option_letter]} km), "
                f"but the provided answer was {llm_answer_letter} ({llm_answer_value} km).")

# Run the check and print the result
result = check_answer()
print(result)