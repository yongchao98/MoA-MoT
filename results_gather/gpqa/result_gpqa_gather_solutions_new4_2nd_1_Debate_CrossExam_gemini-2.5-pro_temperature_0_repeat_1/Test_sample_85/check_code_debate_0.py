def check_chemistry_answer():
    """
    This function checks the correctness of the provided answer to the chemistry problem
    by simulating the chemical logic involved.
    """

    # --- Problem Constraints & Provided Answer ---
    # Reaction A: A + LiBH4 -> (R)-product
    # Reaction B: B + BH3 -> (S)-product
    # Options: A=(S,R), B=(R,S), C=(R,R), D=(S,S)
    provided_answer = 'D'
    product_A_config = 'R'
    product_B_config = 'S'

    # --- Step 1 & 2: Simulate Chemical Logic ---

    # Principle 1: Chemoselectivity
    # LiBH4 reduces the ester. BH3 reduces the acid. This is correct.

    # Principle 2: Stereochemical Consequence based on CIP priority changes.
    # In the starting material, the priority is: Ester side > Acid side.

    # For Reaction A (LiBH4 reduces the ester):
    # The higher-priority ester side becomes a lower-priority alcohol side.
    # The lower-priority acid side remains and becomes the new highest-priority side.
    # The priorities of the top two groups swap.
    # A swap in priority of two groups at a chiral center INVERTS the R/S descriptor.
    stereochem_effect_A = "inversion"

    # For Reaction B (BH3 reduces the acid):
    # The lower-priority acid side becomes an even lower-priority alcohol side.
    # The higher-priority ester side remains the highest-priority side.
    # The priority order is maintained.
    # This means the R/S descriptor is RETAINED.
    stereochem_effect_B = "retention"

    # --- Step 3: Deduce the Correct Starting Materials ---

    # For A: To get an (R) product via INVERSION, the starting material A must be (S).
    # (S)-A --inversion--> (R)-product
    required_A_config = 'S'

    # For B: To get an (S) product via RETENTION, the starting material B must be (S).
    # (S)-B --retention--> (S)-product
    required_B_config = 'S'

    # The derived correct configuration is A=(S) and B=(S).

    # --- Step 4: Match with Options and Verify ---
    options = {
        'A': {'A': 'S', 'B': 'R'},
        'B': {'A': 'R', 'B': 'S'},
        'C': {'A': 'R', 'B': 'R'},
        'D': {'A': 'S', 'B': 'S'}
    }

    derived_correct_option = None
    for option_key, configs in options.items():
        if configs['A'] == required_A_config and configs['B'] == required_B_config:
            derived_correct_option = option_key
            break

    if provided_answer == derived_correct_option:
        return "Correct"
    else:
        reasoning = (
            f"Incorrect. The provided answer is '{provided_answer}', but the correct answer is '{derived_correct_option}'.\n\n"
            "Here is the step-by-step reasoning:\n"
            "1. **Reaction A (A + LiBH₄ → (R)-product):**\n"
            "   - LiBH₄ reduces the ester group. This causes the Cahn-Ingold-Prelog (CIP) priorities of the main side chains to swap.\n"
            "   - A swap in priorities results in an **inversion** of the R/S descriptor.\n"
            "   - To get an (R)-product from an inversion, the starting material A must be **(S)**.\n\n"
            "2. **Reaction B (B + BH₃ → (S)-product):**\n"
            "   - BH₃ reduces the carboxylic acid group. This does not change the CIP priority order of the main side chains.\n"
            "   - No change in priority order results in **retention** of the R/S descriptor.\n"
            "   - To get an (S)-product from a retention, the starting material B must be **(S)**.\n\n"
            "Therefore, both A and B must be the (S)-enantiomer, which corresponds to option D."
        )
        return reasoning

# Run the check
result = check_chemistry_answer()
print(result)