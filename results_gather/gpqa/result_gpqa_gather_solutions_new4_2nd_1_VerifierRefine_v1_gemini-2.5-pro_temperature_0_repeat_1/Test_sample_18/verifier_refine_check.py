import re

def check_chemistry_answer():
    """
    This function checks the correctness of the selected answer for a chemistry question.
    It verifies the products of two Michael addition reactions based on established chemical principles.
    """

    # The final answer provided by the LLM.
    # The provided answer is <<<C>>>, so we extract 'C'.
    final_answer_from_llm = "C"

    # --- Define the properties of the products in each multiple-choice option ---
    # We can parse the product names to determine their key structural features.
    # For Product A: Is the new substituent at position 1 or 3?
    # For Product B: Is the final product a 'succinate' or a 'butanoate' derivative?
    options = {
        "A": {
            "A_name": "methyl 3-(2-((2,4-dimethylphenyl)sulfinyl)ethyl)-2-oxocyclohexane-1-carboxylate",
            "B_name": "4-ethyl 1-methyl 2-cyclopentyl-3,3-diethyl-2-phenylsuccinate"
        },
        "B": {
            "A_name": "methyl 1-(2-((2,4-dimethylphenyl)sulfinyl)ethyl)-2-oxocyclohexane-1-carboxylate",
            "B_name": "4-ethyl 1-methyl 2-cyclopentyl-3,3-diethyl-2-phenylsuccinate"
        },
        "C": {
            "A_name": "methyl 1-(2-((2,4-dimethylphenyl)sulfinyl)ethyl)-2-oxocyclohexane-1-carboxylate",
            "B_name": "ethyl 2-ethyl-2-(1-(2-methoxy-2-oxo-1-phenylethyl)cyclopentyl)butanoate"
        },
        "D": {
            "A_name": "methyl 3-(2-((2,4-dimethylphenyl)sulfinyl)ethyl)-2-oxocyclohexane-1-carboxylate",
            "B_name": "ethyl 2-ethyl-2-(1-(2-methoxy-2-oxo-1-phenylethyl)cyclopentyl)butanoate"
        }
    }

    # --- Define the correct chemical principles for each reaction ---

    # Principle for Reaction A:
    # The Michael donor is methyl 2-oxocyclohexane-1-carboxylate.
    # The proton on the carbon between the two carbonyl groups (C1) is the most acidic.
    # Therefore, the base (NaOEt) will form the enolate at C1, and the Michael addition will occur at this position.
    # The correct product must be a 1-substituted derivative.
    correct_A_substitution_position = 1

    # Principle for Reaction B:
    # The Michael donor is the enolate of ethyl 2-ethylbutanoate.
    # The reaction is a conjugate addition to the Michael acceptor.
    # The resulting product retains the butanoate skeleton from the donor.
    # A 'succinate' structure implies a different carbon skeleton (a 1,4-dicarbonyl) which is not the product of this Michael addition (a 1,5-dicarbonyl).
    correct_B_parent_chain = "butanoate"

    # --- Check the selected answer against the principles ---

    if final_answer_from_llm not in options:
        return f"Invalid Answer: The final answer '{final_answer_from_llm}' is not one of the valid options (A, B, C, D)."

    selected_option = options[final_answer_from_llm]

    # Check Product A's structure
    # We use a simple regex to find the substitution position in the name.
    match_A = re.search(r'methyl (\d)-', selected_option["A_name"])
    if not match_A:
        return f"Could not parse the name for Product A in option {final_answer_from_llm}."
    
    proposed_A_substitution = int(match_A.group(1))

    if proposed_A_substitution != correct_A_substitution_position:
        return (f"Incorrect. The answer '{final_answer_from_llm}' is wrong for Product A. "
                f"Constraint Violated: The Michael addition on methyl 2-oxocyclohexane-1-carboxylate should occur at the C1 position, "
                f"which has the most acidic proton (between two carbonyl groups). "
                f"The selected answer proposes a {proposed_A_substitution}-substituted product, which is incorrect.")

    # Check Product B's structure
    # We check if the name contains 'succinate' or 'butanoate'.
    if "succinate" in selected_option["B_name"].lower():
        proposed_B_type = "succinate"
    elif "butanoate" in selected_option["B_name"].lower():
        proposed_B_type = "butanoate"
    else:
        return f"Could not parse the name for Product B in option {final_answer_from_llm}."

    if proposed_B_type != correct_B_parent_chain:
        return (f"Incorrect. The answer '{final_answer_from_llm}' is wrong for Product B. "
                f"Constraint Violated: The Michael addition of the ethyl 2-ethylbutanoate enolate results in a butanoate derivative. "
                f"The selected answer proposes a '{proposed_B_type}' derivative, which has an incorrect carbon skeleton for this reaction.")

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check and print the result
print(check_chemistry_answer())