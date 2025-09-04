import re

def check_correctness_of_chemistry_answer():
    """
    Checks the correctness of the LLM's answer by applying chemical constraints
    derived from spectroscopic data.
    """
    # 1. Define constraints from the problem description.
    # - FTIR (broad 3000, 1700) + Mass Spec (m/z=45) -> Carboxylic Acid (-COOH)
    # - A carboxylic acid requires 2 oxygen atoms.
    # - A carboxylic acid's C=O bond contributes 1 to the Degree of Unsaturation (DoU).
    required_oxygens = 2
    dou_from_cooh = 1

    # - FTIR (1650) + 1H NMR (vinyl-H) -> Alkene (C=C)
    # - An alkene's C=C bond contributes 1 to the Degree of Unsaturation (DoU).
    dou_from_alkene = 1

    # - Total required DoU is the sum of unsaturations from identified functional groups.
    required_dou = dou_from_cooh + dou_from_alkene

    # 2. Define the candidate formulas and the LLM's answer.
    # The options as listed in the question prompt.
    candidates = {
        "A": "C6H12O",
        "B": "C6H10O2",
        "C": "C6H10O",
        "D": "C6H12O2"
    }
    # The final answer provided by the LLM to be checked.
    llm_answer_option = "B"
    llm_answer_formula = candidates[llm_answer_option]

    # 3. Define helper functions for parsing and calculation.
    def parse_formula(formula_str):
        """Parses a chemical formula string into a dictionary of atom counts."""
        pattern = r'([A-Z])(\d*)'
        matches = re.findall(pattern, formula_str)
        atom_counts = {'C': 0, 'H': 0, 'O': 0}
        for atom, count in matches:
            atom_counts[atom] = int(count) if count else 1
        return atom_counts

    def calculate_dou(counts):
        """Calculates the Degree of Unsaturation for a formula containing C, H, O."""
        C = counts.get('C', 0)
        H = counts.get('H', 0)
        # Formula: DoU = C + 1 - (H/2)
        return C + 1 - (H / 2)

    # 4. Evaluate the LLM's chosen formula against the constraints.
    counts = parse_formula(llm_answer_formula)
    dou = calculate_dou(counts)
    oxygen_count = counts.get('O', 0)

    if oxygen_count != required_oxygens:
        return f"Incorrect. The chosen formula {llm_answer_formula} has {oxygen_count} oxygen(s), but the evidence for a carboxylic acid requires {required_oxygens}."

    if dou != required_dou:
        return f"Incorrect. The chosen formula {llm_answer_formula} has a Degree of Unsaturation (DoU) of {int(dou)}, but the evidence for both a C=O and a C=C bond requires a DoU of {required_dou}."

    # 5. Verify that no other candidate also fits the criteria.
    valid_candidates = []
    for option, formula in candidates.items():
        c = parse_formula(formula)
        if c.get('O', 0) == required_oxygens and calculate_dou(c) == required_dou:
            valid_candidates.append(formula)
    
    if len(valid_candidates) == 1 and valid_candidates[0] == llm_answer_formula:
        return "Correct"
    elif len(valid_candidates) > 1:
        return f"Incorrect. The answer {llm_answer_formula} is a possible solution, but other formulas like {valid_candidates} also satisfy the constraints, making the choice ambiguous."
    else:
        # This case should be caught by the initial checks, but is included for completeness.
        return f"Incorrect. The chosen formula {llm_answer_formula} does not satisfy all constraints."

# Run the check
result = check_correctness_of_chemistry_answer()
print(result)