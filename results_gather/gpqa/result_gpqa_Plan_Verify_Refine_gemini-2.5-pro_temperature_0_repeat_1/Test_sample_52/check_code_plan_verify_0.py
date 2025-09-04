import re

def check_answer():
    """
    Checks the correctness of the LLM's answer by programmatically deriving the
    molecular formula from the spectroscopic data provided in the question.
    """
    # The answer provided by the LLM
    llm_answer_option = "A"
    options = {
        "A": "C11H12O2",
        "B": "C11H14O2",
        "C": "C12H12O2",
        "D": "C12H14O2",
    }
    llm_answer_formula = options.get(llm_answer_option)

    if not llm_answer_formula:
        return f"Invalid answer option '{llm_answer_option}'. Please choose from {list(options.keys())}."

    # --- Step 1: Define molecular fragments based on spectroscopic data ---

    # "di-substituted 6-membered aromatic ring" with "two signals corresponding to aromatic-H"
    # This strongly implies a para-substituted (1,4) benzene ring.
    # Fragment: C6H4
    aromatic_core = {'C': 6, 'H': 4, 'O': 0, 'CH3_count': 0, 'CH2_count': 0}

    # "two signals corresponding to vinyl-H (one doublet and one doublet of quartets)"
    # This is the classic pattern for a propenyl group (-CH=CH-CH3).
    # Fragment: C3H5, which contains one of the two required -CH3 groups.
    propenyl_group = {'C': 3, 'H': 5, 'O': 0, 'CH3_count': 1, 'CH2_count': 0}

    # "FTIR ... shows ... an ester group" and "two signals corresponding to –CH3 groups"
    # and "no signals corresponding to –CH2 groups".
    # This means the second substituent is an ester that contains the second -CH3 group
    # but has no -CH2 groups.
    # Let's define possible ester substituents and check them against the constraints.
    possible_ester_substituents = [
        # e.g., methyl ester group: -C(=O)OCH3
        {'name': 'methyl ester', 'formula': {'C': 2, 'H': 3, 'O': 2}, 'CH3_count': 1, 'CH2_count': 0},
        # e.g., acetate group: -O-C(=O)CH3
        {'name': 'acetate', 'formula': {'C': 2, 'H': 3, 'O': 2}, 'CH3_count': 1, 'CH2_count': 0},
        # An invalid example to test the filter: ethyl ester group -C(=O)OCH2CH3
        {'name': 'ethyl ester', 'formula': {'C': 3, 'H': 5, 'O': 2}, 'CH3_count': 1, 'CH2_count': 1},
    ]

    # --- Step 2: Assemble the molecule and check all constraints ---

    valid_formulas = set()
    substituent1 = propenyl_group

    for substituent2 in possible_ester_substituents:
        # Constraint check: The entire molecule must have no -CH2 groups.
        total_ch2 = aromatic_core['CH2_count'] + substituent1['CH2_count'] + substituent2['CH2_count']
        if total_ch2 > 0:
            continue  # This structure is invalid.

        # Constraint check: The entire molecule must have exactly two -CH3 groups.
        total_ch3 = aromatic_core['CH3_count'] + substituent1['CH3_count'] + substituent2['CH3_count']
        if total_ch3 != 2:
            continue # This structure is invalid.

        # If all constraints are met, calculate the total molecular formula.
        total_C = aromatic_core['C'] + substituent1['C'] + substituent2['formula']['C']
        total_H = aromatic_core['H'] + substituent1['H'] + substituent2['formula']['H']
        total_O = aromatic_core['O'] + substituent1['O'] + substituent2['formula']['O']
        
        formula_string = f"C{total_C}H{total_H}O{total_O}"
        valid_formulas.add(formula_string)

    # --- Step 3: Compare the derived formula(s) with the LLM's answer ---

    if not valid_formulas:
        return "Logic Error: No valid structure could be assembled from the given constraints."

    if llm_answer_formula in valid_formulas:
        return "Correct"
    else:
        expected_formulas = ", ".join(sorted(list(valid_formulas)))
        reason = (
            f"The answer is incorrect. The provided answer corresponds to the formula {llm_answer_formula}. "
            f"However, based on a systematic analysis of the spectroscopic data, the only possible molecular formula is {expected_formulas}.\n"
            f"Reasoning: The structure must contain a C6H4 aromatic core, a C3H5 propenyl group, and an ester fragment (like -COOCH3 or -OCOCH3) with one CH3 and no CH2 groups. "
            f"Summing these fragments (C6H4 + C3H5 + C2H3O2) yields {expected_formulas}."
        )
        return reason

# Run the check
result = check_answer()
print(result)