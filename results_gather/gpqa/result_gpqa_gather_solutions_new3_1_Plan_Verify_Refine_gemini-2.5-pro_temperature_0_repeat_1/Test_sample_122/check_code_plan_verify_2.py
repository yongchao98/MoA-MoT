import math

def check_correctness_of_relativity_answer():
    """
    Checks the correctness of the answer to the special relativity problem.

    The problem involves calculating the distance traveled by a supernova ejecta
    in the Galaxy's reference frame, given the time elapsed in the ejecta's frame.
    This requires using the time dilation formula.
    """

    # --- Problem Parameters ---
    # Velocity of the ejecta (v) in km/s
    v = 60000
    # Proper time in the ejecta's frame (Δt₀) in seconds
    dt0 = 50
    # Speed of light (c) in km/s (using the common approximation)
    c = 300000

    # --- Options from the Question ---
    # A) 2,940,000 km.
    # B) 3,000,000 km.
    # C) 3,060,000 km.
    # D) 2,880,000 km.
    options = {
        "A": 2940000,
        "B": 3000000,
        "C": 3060000,
        "D": 2880000
    }
    
    # The final answer from the LLM to be checked
    llm_answer_key = "C"

    # --- Calculation ---
    # 1. Calculate the Lorentz factor (gamma)
    beta_squared = (v / c)**2
    gamma = 1 / math.sqrt(1 - beta_squared)

    # 2. Calculate the time elapsed in the Galaxy's frame (dilated time)
    dt_galaxy = gamma * dt0

    # 3. Calculate the distance traveled in the Galaxy's frame
    distance_galaxy = v * dt_galaxy

    # --- Verification ---
    # Constraint 1: The relativistic distance must be greater than the classical distance.
    classical_distance = v * dt0
    if options[llm_answer_key] <= classical_distance:
        return (f"Incorrect. The chosen answer {llm_answer_key} ({options[llm_answer_key]:,} km) "
                f"violates a fundamental principle of time dilation. The distance must be "
                f"greater than the classical distance of {classical_distance:,} km (Option B).")

    # Constraint 2: The chosen answer must be the closest to the calculated value.
    # Find the option that is numerically closest to the calculated distance.
    closest_option_key = min(options, key=lambda k: abs(options[k] - distance_galaxy))

    if llm_answer_key == closest_option_key:
        return "Correct"
    else:
        return (f"Incorrect. The calculated distance is approximately {distance_galaxy:,.2f} km. "
                f"The closest option is {closest_option_key} ({options[closest_option_key]:,} km), "
                f"but the provided answer was {llm_answer_key} ({options[llm_answer_key]:,} km).")

# Execute the check and print the result
result = check_correctness_of_relativity_answer()
print(result)