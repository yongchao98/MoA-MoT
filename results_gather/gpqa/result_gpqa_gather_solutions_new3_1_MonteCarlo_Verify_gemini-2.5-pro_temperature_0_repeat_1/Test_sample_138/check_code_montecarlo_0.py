import re

def check_correctness_of_chemistry_answer(final_answer_text):
    """
    Checks the correctness of the final answer for the given chemistry question.

    The core chemical principle is that treating a ketone with (NaNO2, HCl, H2O)
    causes an alpha-oxidation, converting a -CH2- group adjacent to the carbonyl
    into a new carbonyl group, forming a 1,2-diketone.
    Therefore, the starting materials must be ketones.
    """

    # Define the products from the question
    product_A = "4-isopropylcyclohexane-1,2-dione"
    product_B = "5-methylhexane-2,3-dione"

    # Define the options provided in the question
    options = {
        "A": {
            "A": "4-isopropylcyclohexan-1-one",
            "B": "5-methylhexane-2,3-diol"
        },
        "B": {
            "A": "4-isopropylcyclohexan-1-one",
            "B": "5-methylhexan-2-one"
        },
        "C": {
            "A": "4-isopropyl-2-methoxycyclohexan-1-ol",
            "B": "5-methylhexan-2-one"
        },
        "D": {
            "A": "4-isopropyl-2-methoxycyclohexan-1-ol",
            "B": "5-methylhexane-2,3-diol"
        }
    }

    # --- Step 1: Parse the provided answer ---
    match = re.search(r'<<<([A-D])>>>', final_answer_text)
    if not match:
        return f"Invalid answer format. The answer should be in the format <<<X>>> where X is A, B, C, or D. Received: '{final_answer_text}'"
    
    proposed_option_key = match.group(1)
    
    if proposed_option_key not in options:
        return f"Invalid option '{proposed_option_key}'. Please choose from A, B, C, or D."

    proposed_compounds = options[proposed_option_key]
    start_A = proposed_compounds["A"]
    start_B = proposed_compounds["B"]

    # --- Step 2: Determine the correct starting materials based on chemical principles ---
    correct_start_A = "4-isopropylcyclohexan-1-one"
    correct_start_B = "5-methylhexan-2-one"

    # --- Step 3: Check the proposed answer against the correct materials and constraints ---
    
    # Check starting material A from the proposed answer
    a_is_valid = True
    reason_A = ""
    if not start_A.endswith("-one"):
        a_is_valid = False
        reason_A = f"Starting material A ('{start_A}') is incorrect because it is not a ketone (name does not end in '-one'). The reaction requires a ketone substrate."
    elif start_A != correct_start_A:
        a_is_valid = False
        reason_A = f"Starting material A ('{start_A}') is an incorrect ketone. To form '{product_A}', the precursor must be '{correct_start_A}'."

    # Check starting material B from the proposed answer
    b_is_valid = True
    reason_B = ""
    if not start_B.endswith("-one"):
        b_is_valid = False
        reason_B = f"Starting material B ('{start_B}') is incorrect because it is not a ketone (name does not end in '-one'). The reaction requires a ketone substrate."
    elif start_B != correct_start_B:
        b_is_valid = False
        reason_B = f"Starting material B ('{start_B}') is an incorrect ketone. To form '{product_B}', the precursor must be '{correct_start_B}'."

    # --- Step 4: Formulate the final judgment ---
    if a_is_valid and b_is_valid:
        # This means the compounds in the selected option are correct.
        # Now we just need to confirm the option letter is correct.
        if proposed_option_key == "B":
            return "Correct"
        else:
            # This case should not happen if the logic is sound, but as a safeguard:
            return f"The compounds in option '{proposed_option_key}' are correct, but the correct option letter is 'B'."
    else:
        # The selected option is incorrect. Provide the reason.
        error_messages = []
        if not a_is_valid:
            error_messages.append(reason_A)
        if not b_is_valid:
            error_messages.append(reason_B)
        
        return f"Incorrect. The selected option '{proposed_option_key}' is wrong. Reason(s): " + " ".join(error_messages)

# The final answer from the LLM response to be checked
final_answer_from_llm = "<<<B>>>"

# Run the check
result = check_correctness_of_chemistry_answer(final_answer_from_llm)
print(result)