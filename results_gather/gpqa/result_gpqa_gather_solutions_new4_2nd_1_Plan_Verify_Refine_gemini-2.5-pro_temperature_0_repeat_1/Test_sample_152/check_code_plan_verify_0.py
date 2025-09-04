import sys
from io import StringIO

def check_michael_addition_answer():
    """
    This function checks the correctness of the selected answer for the Michael addition question.
    It analyzes each reaction (A, B, and C) based on established chemical principles
    and compares the derived correct compounds with the compounds listed in the chosen option.
    """
    
    # The final answer from the LLM analysis to be checked.
    # The format <<<B>>> is parsed to get the letter 'B'.
    provided_answer_str = "<<<B>>>"
    llm_answer_choice = provided_answer_str.strip('<>').upper()

    # --- Chemical Analysis ---
    # Step 1: Determine the correct product/reactant for each reaction.

    # Reaction A: dimethyl malonate + methyl (E)-3-(p-tolyl)acrylate -> (A)
    # The nucleophile is the enolate of dimethyl malonate, which attacks the beta-carbon of the acrylate.
    # The p-tolyl group is on the beta-carbon, which becomes C2 of the resulting propane backbone.
    # The malonate carbon is C1, and the acrylate alpha-carbon is C3.
    # Structure: (MeOOC)2CH-CH(p-tolyl)-CH2-COOMe
    correct_A = "trimethyl 2-(p-tolyl)propane-1,1,3-tricarboxylate"

    # Reaction B: 1-(cyclohex-1-en-1-yl)piperidine + (E)-but-2-enenitrile -> (B)
    # This is a Stork enamine synthesis. The enamine attacks the beta-carbon of the nitrile.
    # The acidic workup (H3O+) hydrolyzes the intermediate iminium salt back to a ketone.
    # The thermodynamically stable keto form is the major product, not the enol form.
    # The principal functional group is the nitrile, so the parent chain is butanenitrile.
    # The substituent is a (2-oxocyclohexyl) group.
    correct_B = "3-(2-oxocyclohexyl)butanenitrile"

    # Reaction C: C + but-3-en-2-one -> 2-(3-oxobutyl)cyclohexane-1,3-dione
    # This is a retrosynthesis problem. The (3-oxobutyl) group comes from but-3-en-2-one (MVK).
    # It is attached to C2 of the dione ring. This means C2 was the nucleophile.
    # The protons at C2 of a 1,3-dione are highly acidic and are deprotonated by base (KOH) to form the enolate.
    # Therefore, reactant C is the diketone itself.
    correct_C = "cyclohexane-1,3-dione"

    correct_set = {
        "A": correct_A,
        "B": correct_B,
        "C": correct_C
    }

    # Step 2: Define the compound names for each multiple-choice option.
    options = {
        "A": {
            "A": "trimethyl 3-(p-tolyl)propane-1,1,2-tricarboxylate",
            "B": "3-(2-hydroxycyclohex-1-en-1-yl)butanenitrile",
            "C": "cyclohexane-1,3-dione"
        },
        "B": {
            "A": "trimethyl 2-(p-tolyl)propane-1,1,3-tricarboxylate",
            "B": "3-(2-oxocyclohexyl)butanenitrile",
            "C": "cyclohexane-1,3-dione"
        },
        "C": {
            "A": "trimethyl 2-(p-tolyl)propane-1,1,3-tricarboxylate",
            "B": "3-(2-oxocyclohexyl)butanenitrile",
            "C": "2-hydroxycyclohexane-1,3-dione"
        },
        "D": {
            "A": "trimethyl 3-(p-tolyl)propane-1,1,2-tricarboxylate",
            "B": "3-(2-hydroxycyclohex-1-en-1-yl)butanenitrile",
            "C": "2-hydroxycyclohexane-1,3-dione"
        }
    }

    # Step 3: Validate the chosen answer against the correct set.
    if llm_answer_choice not in options:
        return f"Invalid answer choice '{llm_answer_choice}'. The choice must be one of {list(options.keys())}."

    chosen_option_set = options[llm_answer_choice]
    
    errors = []
    # Check component A
    if chosen_option_set["A"] != correct_set["A"]:
        errors.append(f"Component A is incorrect. The answer states A is '{chosen_option_set['A']}', but the correct product is '{correct_set['A']}'. The p-tolyl group should be on the C2 position of the propane backbone, not C3.")
    
    # Check component B
    if chosen_option_set["B"] != correct_set["B"]:
        errors.append(f"Component B is incorrect. The answer states B is '{chosen_option_set['B']}', but the correct product is '{correct_set['B']}'. The acidic workup in a Stork enamine synthesis yields the thermodynamically stable keto tautomer, not the enol form.")

    # Check component C
    if chosen_option_set["C"] != correct_set["C"]:
        errors.append(f"Component C is incorrect. The answer states C is '{chosen_option_set['C']}', but the correct reactant is '{correct_set['C']}'. The Michael donor is the diketone, which forms an enolate at the highly acidic C2 position.")

    if not errors:
        return "Correct"
    else:
        return "Incorrect. Reason(s):\n" + "\n".join(errors)

# Execute the check and print the result
result = check_michael_addition_answer()
print(result)