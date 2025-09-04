import math

def check_answer():
    """
    This function checks the correctness of the provided answer to the binary star system problem.
    """
    # Given data from the question
    # System 1
    P1 = 2.0  # years
    K1_a = 10.0  # km/s
    K1_b = 5.0  # km/s

    # System 2
    P2 = 1.0  # year
    K2_a = 15.0  # km/s
    K2_b = 10.0  # km/s

    # The final answer provided by the LLM
    llm_answer_letter = "C"

    # The options provided in the question
    options = {
        "A": 0.6,
        "B": 0.7,
        "C": 0.4,
        "D": 1.2
    }

    # --- Step 1: Verify the physical principle and formula ---
    # The total mass of a double-lined, eclipsing binary system is proportional to P * (K1 + K2)^3.
    # The ratio of masses M_sys1 / M_sys2 is therefore [P1 * (K1_total)^3] / [P2 * (K2_total)^3].
    # This principle is correctly identified in the provided answer.

    # --- Step 2: Perform the calculation ---
    K_total1 = K1_a + K1_b
    K_total2 = K2_a + K2_b

    # Check if the intermediate sums are correct
    if K_total1 != 15.0:
        return f"Incorrect intermediate calculation: Sum of velocities for system 1 should be 15.0 km/s, but was calculated as {K_total1}."
    if K_total2 != 25.0:
        return f"Incorrect intermediate calculation: Sum of velocities for system 2 should be 25.0 km/s, but was calculated as {K_total2}."

    # Calculate the final ratio
    try:
        calculated_ratio = (P1 * (K_total1)**3) / (P2 * (K_total2)**3)
    except Exception as e:
        return f"An error occurred during calculation: {e}"

    expected_ratio = 54.0 / 125.0  # 0.432
    if not math.isclose(calculated_ratio, expected_ratio, rel_tol=1e-9):
        return f"Calculation error. The calculated ratio is {calculated_ratio:.4f}, but it should be {expected_ratio:.4f} (54/125)."

    # --- Step 3: Compare the result with the options ---
    # Find the option closest to the calculated ratio
    closest_option = min(options.keys(), key=lambda k: abs(options[k] - calculated_ratio))

    # --- Step 4: Check if the LLM's answer is correct ---
    if closest_option == llm_answer_letter:
        return "Correct"
    else:
        return (f"Incorrect final answer. The calculated ratio is {calculated_ratio:.3f}. "
                f"This value is closest to option {closest_option} (~{options[closest_option]}). "
                f"The provided answer was {llm_answer_letter} (~{options[llm_answer_letter]}).")

# Run the check
result = check_answer()
print(result)