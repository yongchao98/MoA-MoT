import re

def check_chemistry_answer(llm_answer):
    """
    Checks the correctness of the answer to the stereochemistry question.

    The function simulates the chemical logic:
    1. Analyzes the selectivity of LiBH4 and BH3.
    2. Analyzes the change in Cahn-Ingold-Prelog (CIP) priorities for each reaction.
    3. Determines if the R/S label inverts or is retained.
    4. Deduces the correct starting materials and compares with the given answer.
    """

    # --- Step 1: Define the problem constraints and options ---
    product_A_config = 'R'
    product_B_config = 'S'
    options = {
        'A': {'A': 'R', 'B': 'R'},
        'B': {'A': 'S', 'B': 'R'},
        'C': {'A': 'S', 'B': 'S'},
        'D': {'A': 'R', 'B': 'S'}
    }

    # --- Step 2: Define groups and their CIP priorities ---
    # We represent the groups and their priorities. A lower index means higher priority.
    # Based on CIP rules: -CH2COOiBu > -CH2COOH > -CH2CH3 > -H
    priority_start = ["ester_side", "acid_side", "ethyl_group"]

    # --- Step 3: Analyze Reaction A (LiBH4) ---
    # LiBH4 reduces the ester to an alcohol.
    # The new groups to compare are the acid and the new alcohol (from the ester).
    # CIP rules: -CH2COOH > -CH2CH2OH.
    # The original #2 priority group (acid_side) is now #1.
    # The original #1 priority group (ester_side) is now #2 (as alcohol_side).
    priority_A_intermediate = ["acid_side", "ester_side_as_alcohol", "ethyl_group"]
    
    # Check if the R/S label inverts. It does because the top two priorities swapped.
    inverts_A = priority_start[0] != priority_A_intermediate[0]

    # Determine the required starting configuration for A
    if inverts_A:
        # To get an (R) product with inversion, we must start with (S).
        required_A_config = 'S' if product_A_config == 'R' else 'R'
    else:
        required_A_config = product_A_config

    # --- Step 4: Analyze Reaction B (BH3) ---
    # BH3 reduces the acid to an alcohol.
    # The new groups to compare are the ester and the new alcohol (from the acid).
    # CIP rules: -CH2COOiBu > -CH2CH2OH.
    # The original #1 priority group (ester_side) is still #1.
    priority_B_intermediate = ["ester_side", "acid_side_as_alcohol", "ethyl_group"]

    # Check if the R/S label inverts. It does not because the top priority group is the same.
    inverts_B = priority_start[0] != priority_B_intermediate[0]

    # Determine the required starting configuration for B
    if inverts_B:
        required_B_config = 'S' if product_B_config == 'S' else 'R'
    else:
        # To get an (S) product with retention, we must start with (S).
        required_B_config = product_B_config

    # --- Step 5: Find the correct option letter ---
    correct_config = {'A': required_A_config, 'B': required_B_config}
    correct_option = None
    for option, configs in options.items():
        if configs == correct_config:
            correct_option = option
            break

    # --- Step 6: Validate the LLM's answer ---
    match = re.search(r'<<<([A-D])>>>', llm_answer)
    if not match:
        return "The answer is in an invalid format. Expected format like <<<C>>>."

    submitted_option = match.group(1)

    if submitted_option == correct_option:
        return "Correct"
    else:
        reason = (
            f"Incorrect. The provided answer is {submitted_option}, but the correct answer is {correct_option}.\n\n"
            f"**Reasoning:**\n"
            f"The key to this problem is understanding how the reaction changes the Cahn-Ingold-Prelog (CIP) priorities of the groups on the chiral center, which determines the R/S label.\n\n"
            f"1.  **Initial Priorities:** In the starting material, the priority order of the main chains is: `-CH2COOiBu` (ester) > `-CH2COOH` (acid).\n\n"
            f"2.  **Reaction A (A + LiBH₄ → (R)-product):**\n"
            f"    - **Constraint (Selectivity):** `LiBH₄` reduces the ester group to an alcohol (`-CH₂CH₂OH`).\n"
            f"    - **Constraint (Stereochemistry):** The priority of the groups changes. The acid group (`-CH₂COOH`) now has higher priority than the new alcohol group. The top two priorities have swapped.\n"
            f"    - **Conclusion:** A swap in priority **inverts** the R/S label. To get an (R)-product, starting material A must be **(S)**.\n\n"
            f"3.  **Reaction B (B + BH₃ → (S)-product):**\n"
            f"    - **Constraint (Selectivity):** `BH₃` reduces the acid group to an alcohol (`-CH₂CH₂OH`).\n"
            f"    - **Constraint (Stereochemistry):** The priority order does not change. The ester group remains the highest priority group.\n"
            f"    - **Conclusion:** Since priorities are maintained, the R/S label is **retained**. To get an (S)-product, starting material B must be **(S)**.\n\n"
            f"**Final correct answer:** A must be (S) and B must be (S), which corresponds to option {correct_option}."
        )
        return reason

# The final answer provided by the LLM to be checked
final_answer_from_llm = "<<<C>>>"

# Run the check
result = check_chemistry_answer(final_answer_from_llm)
print(result)