import math
import re

def check_correctness_of_supernova_answer():
    """
    This function checks the correctness of the final answer for the supernova problem.
    It performs the physics calculation and compares the result to the chosen option.
    """
    # --- Problem Constants and Given Values ---
    v = 60000  # Relative velocity in km/s
    dt0 = 50   # Proper time in seconds
    c = 300000 # Approximate speed of light in km/s

    # --- Options from the Question ---
    options = {
        "A": 2940000,
        "B": 3000000,
        "C": 3060000,
        "D": 2880000
    }

    # The final answer provided by the LLM to be checked
    final_answer_text = "<<<C>>>"

    # --- Step 1: Perform the correct physics calculation ---
    try:
        # Calculate the Lorentz factor (gamma)
        beta_squared = (v / c) ** 2
        if beta_squared >= 1:
            return "Calculation Error: Velocity cannot be >= speed of light."
        gamma = 1 / math.sqrt(1 - beta_squared)

        # Calculate the time elapsed in the Galaxy's frame (dilated time)
        dt_galaxy = gamma * dt0

        # Calculate the distance in the Galaxy's frame
        calculated_distance = v * dt_galaxy
    except Exception as e:
        return f"An error occurred during calculation: {e}"

    # --- Step 2: Check the logical constraint (time dilation) ---
    classical_distance = v * dt0
    
    # --- Step 3: Find the closest option to the calculated result ---
    closest_option_letter = min(options, key=lambda k: abs(options[k] - calculated_distance))

    # --- Step 4: Parse the provided final answer ---
    match = re.search(r'<<<([A-D])>>>', final_answer_text)
    if not match:
        return f"Invalid final answer format: {final_answer_text}"
    chosen_letter = match.group(1)
    chosen_value = options.get(chosen_letter)

    # --- Step 5: Verify the answer ---
    # Check 1: The chosen answer must be greater than the classical distance
    if chosen_value <= classical_distance:
        return (f"Incorrect. The chosen answer {chosen_letter} ({chosen_value} km) violates the principle of time dilation. "
                f"The relativistic distance must be greater than the classical distance of {classical_distance} km.")

    # Check 2: The chosen answer must be the closest to the calculated value
    if chosen_letter != closest_option_letter:
        return (f"Incorrect. The calculated distance is approximately {calculated_distance:.0f} km. "
                f"The closest option is {closest_option_letter} ({options[closest_option_letter]} km), "
                f"but the provided answer was {chosen_letter} ({chosen_value} km).")

    return "Correct"

# The code block returns the result of the check.
# To see the output, you would call the function:
# print(check_correctness_of_supernova_answer())
result = check_correctness_of_supernova_answer()
if result == "Correct":
    print("Correct")
else:
    print(result)