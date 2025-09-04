def check_organic_chemistry_answer():
    """
    This function checks the correctness of the selected answer for a series of Michael addition reactions.
    It codifies the correct products and reactants based on established chemical principles and compares
    them against the chosen option.
    """

    # The final answer provided by the LLM to be checked.
    llm_final_answer = "C"

    # --- Ground Truth Definition ---
    # Based on chemical principles, we define the correct names for A, B, and C.
    # Reaction A: Michael addition of dimethyl malonate to methyl (E)-3-(p-tolyl)acrylate.
    # The nucleophile attacks the β-carbon, leading to a specific isomer.
    correct_product_A = "trimethyl 2-(p-tolyl)propane-1,1,3-tricarboxylate"

    # Reaction B: Stork enamine synthesis followed by acidic workup.
    # The workup hydrolyzes the intermediate to a ketone. The keto form is the major product over the enol form.
    correct_product_B = "3-(2-oxocyclohexyl)butanenitrile"

    # Reaction C: Retro-Michael analysis.
    # The product is 2-(3-oxobutyl)cyclohexane-1,3-dione from but-3-en-2-one.
    # The Michael donor must be the enolate of cyclohexane-1,3-dione.
    correct_reactant_C = "cyclohexane-1,3-dione"

    correct_outcomes = {
        "A": correct_product_A,
        "B": correct_product_B,
        "C": correct_reactant_C
    }

    # --- Problem Space Definition ---
    # The multiple-choice options as provided in the question.
    options = {
        "A": {
            "A": "trimethyl 2-(p-tolyl)propane-1,1,3-tricarboxylate",
            "B": "3-(2-oxocyclohexyl)butanenitrile",
            "C": "2-hydroxycyclohexane-1,3-dione"
        },
        "B": {
            "A": "trimethyl 3-(p-tolyl)propane-1,1,2-tricarboxylate",
            "B": "3-(2-hydroxycyclohex-1-en-1-yl)butanenitrile",
            "C": "cyclohexane-1,3-dione"
        },
        "C": {
            "A": "trimethyl 2-(p-tolyl)propane-1,1,3-tricarboxylate",
            "B": "3-(2-oxocyclohexyl)butanenitrile",
            "C": "cyclohexane-1,3-dione"
        },
        "D": {
            "A": "trimethyl 3-(p-tolyl)propane-1,1,2-tricarboxylate",
            "B": "3-(2-hydroxycyclohex-1-en-1-yl)butanenitrile",
            "C": "2-hydroxycyclohexane-1,3-dione"
        }
    }

    # --- Verification Logic ---
    if llm_final_answer not in options:
        return f"Invalid answer choice '{llm_final_answer}'. The valid options are A, B, C, D."

    chosen_option_details = options[llm_final_answer]
    
    # Check each part of the chosen answer against the ground truth.
    if chosen_option_details["A"] != correct_outcomes["A"]:
        return (f"Incorrect. The answer '{llm_final_answer}' is wrong because of component A. "
                f"It states product A is '{chosen_option_details['A']}', but the correct product is '{correct_outcomes['A']}'. "
                "The Michael addition places the p-tolyl group on the C2 of the resulting propane-tricarboxylate backbone.")

    if chosen_option_details["B"] != correct_outcomes["B"]:
        return (f"Incorrect. The answer '{llm_final_answer}' is wrong because of component B. "
                f"It states product B is '{chosen_option_details['B']}', but the correct product is '{correct_outcomes['B']}'. "
                "The acidic workup in the Stork enamine synthesis yields the more stable keto form ('oxo'), not the enol form ('hydroxy...en').")

    if chosen_option_details["C"] != correct_outcomes["C"]:
        return (f"Incorrect. The answer '{llm_final_answer}' is wrong because of component C. "
                f"It states reactant C is '{chosen_option_details['C']}', but the correct reactant is '{correct_outcomes['C']}'. "
                "The Michael donor that forms the product is cyclohexane-1,3-dione, which is deprotonated at its acidic C2 position.")

    # If all components match the correct outcomes, the answer is correct.
    return "Correct"

# Execute the check and print the result.
result = check_organic_chemistry_answer()
print(result)