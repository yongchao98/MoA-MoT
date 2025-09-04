import re

def check_correctness():
    """
    Checks the correctness of the LLM's answer by verifying the chemical formula against
    constraints derived from the provided spectroscopic data.
    """
    # --- 1. Define Constraints from the Question ---
    # Evidence for -COOH group implies 2 oxygen atoms.
    required_oxygens = 2
    # Evidence for -COOH (C=O) and an alkene (C=C) implies a total of 2 degrees of unsaturation.
    required_dou = 2

    # --- 2. Define the Options and the LLM's Answer ---
    options = {
        "A": "C6H10O",
        "B": "C6H10O2",
        "C": "C6H12O2",
        "D": "C6H12O"
    }
    llm_answer = "B" # Extracted from <<<B>>>

    # --- 3. Helper Functions ---
    def parse_formula(formula_str):
        """Parses a chemical formula string into a dictionary of element counts."""
        counts = {'C': 0, 'H': 0, 'O': 0}
        # Find all element-count pairs (e.g., C6, H10, O2)
        for element, count in re.findall(r'([A-Z])(\d*)', formula_str):
            counts[element] = int(count) if count else 1
        return counts

    def calculate_dou(C, H):
        """Calculates the Degree of Unsaturation for a C, H, O compound."""
        return C + 1 - (H / 2)

    # --- 4. Evaluate Each Option ---
    valid_options = []
    reasons_for_failure = {}

    for option_letter, formula in options.items():
        counts = parse_formula(formula)
        C, H, O = counts['C'], counts['H'], counts['O']
        
        dou = calculate_dou(C, H)

        # Check if all constraints are met
        is_oxygen_correct = (O == required_oxygens)
        is_dou_correct = (dou == required_dou)

        if is_oxygen_correct and is_dou_correct:
            valid_options.append(option_letter)
        else:
            # Record why an option failed for better error messages
            failures = []
            if not is_oxygen_correct:
                failures.append(f"it has {O} oxygen(s), but the evidence for a carboxylic acid requires {required_oxygens}")
            if not is_dou_correct:
                failures.append(f"its Degree of Unsaturation is {int(dou)}, but the evidence for a C=O and a C=C bond requires a DoU of {required_dou}")
            reasons_for_failure[option_letter] = " and ".join(failures)

    # --- 5. Final Verdict ---
    if not valid_options:
        return "Incorrect. The reasoning is flawed as none of the options satisfy all the spectral data constraints."
    
    if len(valid_options) > 1:
        return f"Incorrect. The problem is ambiguous as options {', '.join(valid_options)} both satisfy the constraints."

    correct_option = valid_options[0]

    if llm_answer == correct_option:
        return "Correct"
    else:
        reason = reasons_for_failure.get(llm_answer, "it does not meet the derived constraints")
        return f"Incorrect. The provided answer is {llm_answer}, but the only option that satisfies all constraints is {correct_option}. The answer {llm_answer} is wrong because {reason}."

# Execute the check and print the result
result = check_correctness()
print(result)