import re

def check_answer():
    """
    Checks the correctness of the selected chemical formula by systematically
    analyzing the constraints from the spectroscopic data.
    """
    # --- Problem Data & LLM's Final Answer ---
    options = {
        'A': 'C12H14O2',
        'B': 'C11H12O2',
        'C': 'C11H14O2',
        'D': 'C12H12O2'
    }
    # The final answer from the provided text is 'B'.
    llm_answer_key = 'B'
    llm_answer_formula = options.get(llm_answer_key)

    if not llm_answer_formula:
        return f"Invalid answer key '{llm_answer_key}'. Please use one of {list(options.keys())}."

    # --- Step 1: Deduce the correct formula from the spectroscopic constraints ---

    # Constraint: Di-substituted 6-membered aromatic ring -> C6H4 fragment
    # Constraint: 1H NMR shows two vinyl-H (doublet, doublet of quartets) -> propenyl group (-CH=CH-CH3) -> C3H5 fragment
    # Constraint: FTIR shows ester group (-COO-), 1H NMR shows two -CH3 signals, and no -CH2 groups.
    # This combination points to a methyl ester (-COOCH3) or acetate (-OCOCH3). Both add C2H3O2 to the formula.
    
    fragments = {
        'aromatic_ring': {'C': 6, 'H': 4, 'O': 0},
        'propenyl_group': {'C': 3, 'H': 5, 'O': 0},
        'ester_with_methyl': {'C': 2, 'H': 3, 'O': 2}
    }

    # Sum the atoms from all fragments to get the expected formula
    expected_counts = {'C': 0, 'H': 0, 'O': 0}
    for frag_counts in fragments.values():
        for elem, count in frag_counts.items():
            expected_counts[elem] += count
    
    expected_formula = f"C{expected_counts['C']}H{expected_counts['H']}O{expected_counts['O']}"

    # --- Step 2: Verify the LLM's answer against the deduced formula ---
    if llm_answer_formula != expected_formula:
        return (f"Incorrect. The provided answer is {llm_answer_key} ({llm_answer_formula}), but the formula deduced from the spectral data is {expected_formula}. "
                f"The fragments (C6H4 ring, C3H5 propenyl, and C2H3O2 ester part) sum to {expected_formula}.")

    # --- Step 3: Verify using Degree of Unsaturation (DoU) as a cross-check ---
    def parse_formula(formula_str):
        counts = re.findall(r'([A-Z])(\d*)', formula_str)
        return {elem: int(count) if count else 1 for elem, count in counts}

    def calculate_dou(counts):
        c = counts.get('C', 0)
        h = counts.get('H', 0)
        return c + 1 - (h / 2)

    # Required DoU from structure: 4 (aromatic ring) + 1 (C=O ester) + 1 (C=C vinyl) = 6
    required_dou = 6.0
    
    llm_answer_counts = parse_formula(llm_answer_formula)
    if calculate_dou(llm_answer_counts) != required_dou:
        return (f"Incorrect. The Degree of Unsaturation for the answer {llm_answer_formula} is {calculate_dou(llm_answer_counts)}, "
                f"but the required DoU based on the structure (aromatic ring + C=O + C=C) is {required_dou}.")

    # --- Step 4: Explain why other options are incorrect ---
    # Option A: C12H14O2. DoU = 12 + 1 - (14/2) = 6. Same DoU as the correct answer.
    # However, C12H14O2 is equivalent to C11H12O2 + CH2. This addition of a CH2 group
    # violates the "no signals corresponding to –CH2 groups" constraint.
    
    # Option C: C11H14O2. DoU = 11 + 1 - (14/2) = 5. This DoU is incorrect.
    
    # Option D: C12H12O2. DoU = 12 + 1 - (12/2) = 7. This DoU is incorrect.

    # The LLM's answer is the only one that satisfies all constraints.
    return "Correct"

# Run the check
result = check_answer()
print(result)