import math

def check_answer():
    """
    Checks the correctness of the answer to the special relativity problem.
    """
    # --- Problem Parameters ---
    # Velocity of the ejecta relative to the Galaxy
    v = 60000  # km/s
    # Proper time (time elapsed in the ejecta's frame)
    t0 = 50    # s
    # Speed of light (approximated as is standard for such problems)
    c = 300000 # km/s

    # --- Options from the Question ---
    # A) 2 940 000 km.
    # B) 3 060 000 km.
    # C) 3 000 000 km.
    # D) 2 880 000 km.
    options = {
        "A": 2940000,
        "B": 3060000,
        "C": 3000000,
        "D": 2880000
    }
    
    # The final answer provided by the LLM to be checked
    proposed_answer_letter = "B"

    # --- Relativistic Calculation ---
    # Step 1: Calculate the Lorentz factor (gamma)
    # gamma = 1 / sqrt(1 - v^2/c^2)
    try:
        beta_squared = (v / c)**2
        if beta_squared >= 1:
            return "Calculation error: Velocity cannot be equal to or greater than the speed of light."
        gamma = 1 / math.sqrt(1 - beta_squared)
    except Exception as e:
        return f"An error occurred during calculation: {e}"

    # Step 2: Calculate the time elapsed in the Galaxy's frame (dilated time)
    # t_galaxy = gamma * t0
    t_galaxy = gamma * t0

    # Step 3: Calculate the distance traveled in the Galaxy's frame
    # d = v * t_galaxy
    calculated_distance = v * t_galaxy

    # --- Verification ---
    # Find which option is numerically closest to the calculated distance
    closest_option_letter = None
    min_difference = float('inf')

    for letter, value in options.items():
        difference = abs(calculated_distance - value)
        if difference < min_difference:
            min_difference = difference
            closest_option_letter = letter

    # Check if the proposed answer matches the closest option
    if proposed_answer_letter == closest_option_letter:
        return "Correct"
    else:
        # Provide a detailed reason for the incorrectness
        reason = (
            f"The proposed answer '{proposed_answer_letter}' ({options[proposed_answer_letter]:,} km) is incorrect.\n"
            f"The correct calculation based on special relativity yields a distance of {calculated_distance:,.2f} km.\n\n"
            f"Calculation Steps:\n"
            f"1. Lorentz Factor (γ) = 1 / sqrt(1 - (v/c)²) = {gamma:.6f}\n"
            f"2. Dilated Time in Galaxy Frame (t) = γ * t₀ = {t_galaxy:.4f} s\n"
            f"3. Distance in Galaxy Frame (d) = v * t = {calculated_distance:,.2f} km\n\n"
            f"This calculated distance is closest to option '{closest_option_letter}' ({options[closest_option_letter]:,} km).\n"
            f"A common mistake is to ignore relativity and calculate d = v * t₀ = {v*t0:,} km, which corresponds to option 'C'."
        )
        return reason

# Run the check and print the result
print(check_answer())