def check_diels_alder_product():
    """
    This function checks the correctness of the proposed answer by verifying
    if the product's stereochemistry satisfies the constraints of the
    Diels-Alder reaction described in the question.
    """

    # --- Step 1: Define the constraints from the reaction ---

    # Constraint 1: The product must have a 4,5,6-trimethylcyclohex-2-enol skeleton.
    # This comes from the Diels-Alder reaction of penta-1,3-dien-1-ol and but-2-ene.
    constraint_skeleton_is_456 = True

    # Constraint 2: The methyl groups at C5 and C6 must be cis to each other.
    # This is because the dienophile is cis-but-2-ene.
    constraint_c5_c6_must_be_cis = True

    # Constraint 3: The hydroxyl group at C1 and the methyl group at C4 must be cis to each other.
    # This is because the diene is (1E,3E)-penta-1,3-dien-1-ol, where the terminal
    # substituents are on the "outside" of the s-cis conformer.
    constraint_c1_c4_must_be_cis = True

    # --- Step 2: Analyze the stereochemistry of the given options ---
    # This analysis is done manually by converting R/S notation to relative
    # cis/trans relationships for the specific molecular structure.
    # Let's assume 'down' and 'up' for substituent orientation.
    # For (1S, ...), we can set the C1-OH group to be 'down'.

    options = {
        "A": {
            "name": "(1S,4R,5S,6R)-4,5,6-trimethylcyclohex-2-enol",
            "has_correct_skeleton": True,
            # Analysis: C1-OH(down), C4-Me(up), C5-Me(up), C6-Me(up)
            "is_c5_c6_cis": True,  # Both C5-Me and C6-Me are 'up'
            "is_c1_c4_cis": False, # C1-OH is 'down', C4-Me is 'up'
        },
        "B": {
            "name": "(1S,4S)-4,6,6-trimethylcyclohex-2-enol",
            "has_correct_skeleton": False, # 4,6,6-trimethyl is the wrong skeleton
            "is_c5_c6_cis": None, # Not applicable
            "is_c1_c4_cis": None, # Not applicable
        },
        "C": {
            "name": "(1S,4R)-4,6,6-trimethylcyclohex-2-enol",
            "has_correct_skeleton": False, # 4,6,6-trimethyl is the wrong skeleton
            "is_c5_c6_cis": None, # Not applicable
            "is_c1_c4_cis": None, # Not applicable
        },
        "D": {
            "name": "(1S,4R,5S,6S)-4,5,6-trimethylcyclohex-2-enol",
            "has_correct_skeleton": True,
            # Analysis: C1-OH(down), C4-Me(up), C5-Me(up), C6-Me(down)
            "is_c5_c6_cis": False, # C5-Me is 'up', C6-Me is 'down'
            "is_c1_c4_cis": False, # C1-OH is 'down', C4-Me is 'up'
        }
    }

    llm_answer = "A"
    chosen_option = options[llm_answer]

    # --- Step 3: Check the LLM's answer against the constraints ---
    
    error_messages = []

    if not chosen_option["has_correct_skeleton"]:
        error_messages.append(f"The product skeleton is incorrect. Expected a 4,5,6-trimethyl structure, but the answer has a {chosen_option['name']}.")
    
    if chosen_option["is_c5_c6_cis"] != constraint_c5_c6_must_be_cis:
        error_messages.append("Constraint from dienophile (cis-but-2-ene) is not satisfied. The methyl groups at C5 and C6 must be cis, but they are trans in the answer.")

    if chosen_option["is_c1_c4_cis"] != constraint_c1_c4_must_be_cis:
        error_messages.append("Constraint from diene ((1E,3E)-penta-1,3-dien-1-ol) is not satisfied. The -OH group at C1 and the methyl group at C4 must be cis, but they are trans in the answer.")

    if not error_messages:
        return "Correct"
    else:
        return f"The answer '{llm_answer}' is incorrect. Reason(s):\n- " + "\n- ".join(error_messages)

# Run the check
result = check_diels_alder_product()
print(result)