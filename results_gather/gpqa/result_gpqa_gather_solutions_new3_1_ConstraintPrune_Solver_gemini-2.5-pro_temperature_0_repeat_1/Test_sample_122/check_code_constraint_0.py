import math

def check_correctness_of_supernova_answer():
    """
    This function checks the correctness of the provided answer to the special relativity problem.

    The problem asks for the distance traveled by a supernova ejecta in the Galaxy's reference frame.
    - Velocity (v) = 60,000 km/s
    - Proper time (Δt₀) = 50 s
    - Speed of light (c) ≈ 300,000 km/s

    The solution involves calculating the time dilation effect.
    1. Calculate the Lorentz factor (γ) = 1 / sqrt(1 - v²/c²)
    2. Calculate the time in the Galaxy's frame (Δt) = γ * Δt₀
    3. Calculate the distance in the Galaxy's frame (d) = v * Δt
    """

    # --- Problem Parameters ---
    v = 60000  # km/s
    delta_t0 = 50  # seconds (proper time)
    c = 300000  # km/s

    # --- Options from the Question ---
    # A) 2 880 000 km.
    # B) 3 000 000 km.
    # C) 3 060 000 km.
    # D) 2 940 000 km.
    options = {
        "A": 2880000,
        "B": 3000000,
        "C": 3060000,
        "D": 2940000
    }

    # The final answer provided by the LLM
    llm_final_answer_letter = "C"
    
    # --- Physics Calculation ---
    try:
        # Calculate the Lorentz factor (gamma)
        beta_squared = (v / c) ** 2
        gamma = 1 / math.sqrt(1 - beta_squared)

        # Calculate the time elapsed in the Galaxy's frame (dilated time)
        delta_t_galaxy = gamma * delta_t0

        # Calculate the distance traveled as measured in the Galaxy's frame
        calculated_distance = v * delta_t_galaxy
    except (ValueError, ZeroDivisionError) as e:
        return f"An error occurred during calculation: {e}"

    # --- Verification ---
    # Find which option is closest to the calculated distance
    closest_option = None
    min_difference = float('inf')

    for letter, value in options.items():
        difference = abs(calculated_distance - value)
        if difference < min_difference:
            min_difference = difference
            closest_option = letter

    # Check if the LLM's answer matches the calculated closest option
    if llm_final_answer_letter == closest_option:
        return "Correct"
    else:
        reason = (
            f"The provided answer is '{llm_final_answer_letter}', which corresponds to {options[llm_final_answer_letter]:,} km.\n"
            f"However, the correct calculation yields a different result.\n\n"
            f"Step-by-step calculation:\n"
            f"1. Lorentz factor (γ) = 1 / sqrt(1 - (60000/300000)²) = 1 / sqrt(0.96) ≈ {gamma:.6f}\n"
            f"2. Time in Galaxy frame (Δt) = γ * 50s ≈ {delta_t_galaxy:.4f} s\n"
            f"3. Distance = v * Δt = 60,000 km/s * {delta_t_galaxy:.4f} s ≈ {calculated_distance:,.2f} km\n\n"
            f"The calculated distance ({calculated_distance:,.2f} km) is closest to option '{closest_option}' ({options[closest_option]:,} km), not '{llm_final_answer_letter}'."
        )
        return reason

# Execute the check and print the result
print(check_correctness_of_supernova_answer())