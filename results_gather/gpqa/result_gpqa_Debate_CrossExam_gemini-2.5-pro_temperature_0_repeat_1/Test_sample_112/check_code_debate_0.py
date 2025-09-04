import re

def check_correctness():
    """
    Checks the correctness of the LLM's answer by verifying which chemical formula
    matches the constraints derived from the spectroscopic data.
    """
    # The final answer provided by the LLM
    llm_answer_key = "D"

    # The options given in the question
    options = {
        "A": "C6H12O",
        "B": "C6H10O",
        "C": "C6H12O2",
        "D": "C6H10O2"
    }

    # --- Constraints derived from the spectroscopic data ---
    # 1. FTIR (broad 3000, 1700) & Mass Spec (m/z=45) indicate a carboxylic acid (-COOH).
    #    This requires exactly 2 oxygen atoms.
    required_oxygens = 2

    # 2. FTIR (1650) & 1H NMR (vinyl-H) indicate a C=C double bond.
    #    The carboxylic acid has a C=O bond.
    #    Total Degree of Unsaturation (DoU) = 1 (from C=O) + 1 (from C=C) = 2.
    required_dou = 2

    def parse_formula(formula_str):
        """Parses a chemical formula string into a dictionary of atom counts."""
        # Finds all element-count pairs (e.g., C6, H12, O)
        parts = re.findall(r'([A-Z])(\d*)', formula_str)
        counts = {'C': 0, 'H': 0, 'O': 0}
        for element, count in parts:
            counts[element] = int(count) if count else 1
        return counts

    def calculate_dou(c, h):
        """Calculates the Degree of Unsaturation for a CxHyOz formula."""
        # Formula: DoU = C - H/2 + 1
        return c - (h / 2) + 1

    # --- Evaluate each option against the constraints ---
    valid_options = []
    reasons_for_failure = {}

    for key, formula in options.items():
        counts = parse_formula(formula)
        c, h, o = counts['C'], counts['H'], counts['O']

        # Check 1: Number of oxygen atoms
        if o != required_oxygens:
            reasons_for_failure[key] = f"Fails the oxygen count constraint. The evidence for a carboxylic acid requires {required_oxygens} oxygen atoms, but formula {formula} has {o}."
            continue

        # Check 2: Degree of Unsaturation
        dou = calculate_dou(c, h)
        if dou != required_dou:
            reasons_for_failure[key] = f"Fails the Degree of Unsaturation (DoU) constraint. The evidence for a C=C and a C=O group requires a DoU of {required_dou}, but formula {formula} has a DoU of {int(dou)}."
            continue

        # If a formula passes all checks, it's a valid candidate
        valid_options.append(key)

    # --- Final conclusion ---
    if len(valid_options) == 1 and valid_options[0] == llm_answer_key:
        return "Correct"
    elif len(valid_options) != 1:
        return f"The analysis is flawed. Found {len(valid_options)} valid options: {valid_options}. A single correct answer should be identifiable."
    else:
        correct_key = valid_options[0]
        llm_reason_for_failure = reasons_for_failure.get(llm_answer_key, "Unknown failure.")
        return (f"Incorrect. The LLM's answer is {llm_answer_key}, but the only option that satisfies all constraints is {correct_key} ({options[correct_key]}).\n"
                f"Reason why the LLM's choice ({options[llm_answer_key]}) is wrong: {llm_reason_for_failure}")

# Run the check
result = check_correctness()
print(result)