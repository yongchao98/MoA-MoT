import re

def check_correctness():
    """
    Checks the correctness of the LLM's answer by applying chemical constraints
    derived from the spectroscopic data.
    """
    # --- 1. Define constraints from the problem description ---
    # Evidence for Carboxylic Acid (-COOH):
    # - FTIR: Broad peak at 3000 cm-1 (O-H), strong peak at 1700 cm-1 (C=O)
    # - Mass Spec: Fragment at m/z = 45 ([COOH]+)
    # This implies the molecule must have 2 oxygen atoms and 1 DoU from C=O.
    required_oxygens = 2
    dou_from_carbonyl = 1

    # Evidence for Alkene (C=C):
    # - FTIR: Strong peak at 1650 cm-1
    # - 1H NMR: Vinyl-hydrogens
    # This implies the molecule must have 1 DoU from C=C.
    dou_from_alkene = 1

    # Total required Degree of Unsaturation (DoU)
    required_dou = dou_from_carbonyl + dou_from_alkene

    # --- 2. Define the options and the LLM's answer ---
    # The options as presented in the original question text.
    options = {
        "A": "C6H10O2",
        "B": "C6H10O",
        "C": "C6H12O",
        "D": "C6H12O2"
    }
    # The final answer provided by the LLM to be checked.
    llm_answer = "A"

    # --- 3. Helper functions for chemical formula analysis ---
    def get_atom_count(formula, atom_symbol):
        """Parses a chemical formula string to find the count of a specific atom."""
        match = re.search(f"{atom_symbol}(\d*)", formula)
        if not match:
            return 0
        count_str = match.group(1)
        return int(count_str) if count_str else 1

    def calculate_dou(formula):
        """Calculates the Degree of Unsaturation for a C, H, O formula."""
        num_c = get_atom_count(formula, 'C')
        num_h = get_atom_count(formula, 'H')
        # DoU = C + 1 - (H/2)
        return num_c + 1 - (num_h / 2)

    # --- 4. Evaluate all candidates and find the correct one ---
    valid_options = []
    failure_reasons = {}

    for letter, formula in options.items():
        is_valid = True
        reasons = []

        # Check Constraint 1: Oxygen count
        oxygen_count = get_atom_count(formula, 'O')
        if oxygen_count != required_oxygens:
            is_valid = False
            reasons.append(f"has {oxygen_count} oxygen(s), but a carboxylic acid requires {required_oxygens}")

        # Check Constraint 2: Degree of Unsaturation
        dou = calculate_dou(formula)
        if dou != required_dou:
            is_valid = False
            reasons.append(f"has a Degree of Unsaturation (DoU) of {int(dou)}, but the evidence (C=O and C=C) requires a DoU of {required_dou}")

        if is_valid:
            valid_options.append(letter)
        else:
            failure_reasons[letter] = " and ".join(reasons)

    # --- 5. Compare the derived correct answer with the LLM's answer ---
    if len(valid_options) != 1:
        return f"Analysis Error: The problem constraints lead to {len(valid_options)} valid solutions ({', '.join(valid_options)}), not a unique answer."

    correct_option = valid_options[0]

    if llm_answer == correct_option:
        return "Correct"
    else:
        reason_for_llm_failure = failure_reasons.get(llm_answer, "it does not meet the derived constraints")
        return (f"Incorrect. The provided answer is '{llm_answer}', but the correct answer is '{correct_option}'. "
                f"The formula for option '{correct_option}' ({options[correct_option]}) is the only one that satisfies all constraints. "
                f"The chosen answer '{llm_answer}' ({options[llm_answer]}) is wrong because it {reason_for_llm_failure}.")

# Run the check and print the result
result = check_correctness()
print(result)