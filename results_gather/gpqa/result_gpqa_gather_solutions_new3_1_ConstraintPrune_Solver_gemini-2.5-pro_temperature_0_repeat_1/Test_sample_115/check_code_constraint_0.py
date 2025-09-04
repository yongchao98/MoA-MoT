import re

def check_chemistry_answer():
    """
    Checks the correctness of the LLM's answer by applying chemical constraints.
    """
    # --- Problem Definition ---
    question_data = {
        "dienophile": "cis-but-2-ene",
        "diene": "(1E,3E)-penta-1,3-dien-1-ol",
        "options": {
            "A": "(1S,4R,5S,6R)-4,5,6-trimethylcyclohex-2-enol",
            "B": "(1S,4S)-4,6,6-trimethylcyclohex-2-enol",
            "C": "(1S,4R,5S,6S)-4,5,6-trimethylcyclohex-2-enol",
            "D": "(1S,4R)-4,6,6-trimethylcyclohex-2-enol"
        },
        "llm_answer_key": "A"
    }

    # --- Constraint 1: Regiochemistry ---
    # The product must be a 4,5,6-trimethylcyclohex-2-enol.
    valid_options_regio = {}
    for key, value in question_data["options"].items():
        if "4,5,6-trimethyl" in value:
            valid_options_regio[key] = value

    # Check if the LLM's answer satisfies this constraint
    if question_data["llm_answer_key"] not in valid_options_regio:
        return (f"Incorrect: The provided answer '{question_data['llm_answer_key']}' has the wrong regiochemistry. "
                f"The product skeleton should be '4,5,6-trimethylcyclohex-2-enol', but the answer corresponds to a "
                f"'{question_data['options'][question_data['llm_answer_key']]}'.")

    # --- Constraint 2: Dienophile Stereospecificity ---
    # The dienophile is cis, so substituents at C5 and C6 must be cis.
    # For adjacent stereocenters, cis means (R,S) or (S,R). Trans means (R,R) or (S,S).
    final_valid_options = {}
    for key, value in valid_options_regio.items():
        # Use regex to find the stereochemistry at C5 and C6
        match = re.search(r'5(S|R),6(S|R)', value)
        if not match:
            # This option doesn't have stereocenters at C5 and C6, so it's incorrect.
            continue
        
        c5_config = match.group(1)
        c6_config = match.group(2)

        # A cis relationship means the R/S descriptors are different.
        if c5_config != c6_config:
            final_valid_options[key] = value

    # --- Final Verification ---
    if not final_valid_options:
        return ("Incorrect: No option satisfies both the regiochemistry and the dienophile stereospecificity rules.")

    if len(final_valid_options) == 1 and question_data["llm_answer_key"] in final_valid_options:
        return "Correct"
    else:
        # The LLM's answer was eliminated. Explain why.
        if question_data["llm_answer_key"] in valid_options_regio:
            # It passed regiochemistry but failed stereochemistry
            return (f"Incorrect: The provided answer '{question_data['llm_answer_key']}' violates the dienophile stereospecificity rule. "
                    f"Since the dienophile is cis-but-2-ene, the substituents at C5 and C6 must be cis (i.e., have opposite R/S descriptors like 5S,6R). "
                    f"The answer has a trans relationship at C5 and C6.")
        else:
            # This case is already handled by the first check, but included for completeness.
            return f"Incorrect: The provided answer '{question_data['llm_answer_key']}' was eliminated by the constraint checks."


# Run the check
result = check_chemistry_answer()
print(result)