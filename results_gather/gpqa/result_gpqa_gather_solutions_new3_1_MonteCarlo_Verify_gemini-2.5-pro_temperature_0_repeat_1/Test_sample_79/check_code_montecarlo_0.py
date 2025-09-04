def check_answer():
    """
    Checks the logic and calculations of the provided answer.
    """
    # --- Step 1: Define the problem and the proposed solution ---
    options = {'A': 315, 'B': 333, 'C': 185, 'D': 351}
    proposed_answer_letter = 'B'
    proposed_answer_value = options[proposed_answer_letter]

    # The specific values found in the solution to verify the logic
    solution_values = {'A': 25, 'C': 55, 'G': 45, 'T': 3}

    # --- Step 2: Verify the specific values found in the solution ---
    A = solution_values['A']
    C = solution_values['C']
    G = solution_values['G']
    T = solution_values['T']

    # Check if the values satisfy the first equation from "AGG -> 115"
    if A + 2 * G != 115:
        return f"Incorrect: The proposed values {solution_values} do not satisfy the first constraint A + 2G = 115. Got {A + 2 * G}."

    # Check if the values satisfy the second equation from "TGCTGA -> 176"
    # This simplifies to C + 2T = 61
    if C + 2 * T != 61:
        return f"Incorrect: The proposed values {solution_values} do not satisfy the second constraint C + 2T = 61. Got {C + 2 * T}."

    # Check if the values produce the proposed answer for "ACAGTGACC"
    # Target expression: 3A + 3C + 2G + T
    calculated_value = 3 * A + 3 * C + 2 * G + T
    if calculated_value != proposed_answer_value:
        return f"Incorrect: The proposed values {solution_values} calculate to {calculated_value}, not the proposed answer {proposed_answer_value}."

    # --- Step 3: Verify the "hidden constraint" logic ---
    # The logic is that only S=333 forces G to be a multiple of 5.
    # The derived equation is 4G + 5T = 528 - S.
    # If 528 - S is a multiple of 5, then 4G must be, and thus G must be.
    
    g_must_be_multiple_of_5_count = 0
    option_with_constraint = None

    for letter, s_value in options.items():
        k = 528 - s_value
        if k % 5 == 0:
            g_must_be_multiple_of_5_count += 1
            option_with_constraint = letter
    
    if g_must_be_multiple_of_5_count != 1:
        return f"Incorrect: The hidden constraint logic is flawed. Found {g_must_be_multiple_of_5_count} options that force G to be a multiple of 5, not just one."

    if option_with_constraint != proposed_answer_letter:
        return f"Incorrect: The hidden constraint (G is a multiple of 5) points to option {option_with_constraint}, not the proposed answer {proposed_answer_letter}."

    # --- Step 4: Final check ---
    # All checks passed, the logic is sound and the calculations are correct.
    return "Correct"

# Run the checker
result = check_answer()
print(result)