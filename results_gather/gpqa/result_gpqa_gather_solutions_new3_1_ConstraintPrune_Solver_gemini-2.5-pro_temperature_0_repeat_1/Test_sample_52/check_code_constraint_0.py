import re

def parse_formula(formula):
    """Parses a chemical formula string into a dictionary of atom counts."""
    pattern = r'([A-Z][a-z]*)(\d*)'
    atoms = {}
    for element, count in re.findall(pattern, formula):
        atoms[element] = int(count) if count else 1
    return atoms

def calculate_dou(atoms):
    """Calculates the Degree of Unsaturation (Double Bond Equivalent)."""
    c = atoms.get('C', 0)
    h = atoms.get('H', 0)
    n = atoms.get('N', 0)
    x = atoms.get('F', 0) + atoms.get('Cl', 0) + atoms.get('Br', 0) + atoms.get('I', 0)
    # Oxygen and sulfur do not affect the calculation
    return c - (h / 2) - (x / 2) + (n / 2) + 1

def check_answer():
    """
    Checks the correctness of the answer based on the problem's constraints.
    """
    question_options = {
        "A": "C12H12O2",
        "B": "C11H12O2",
        "C": "C12H14O2",
        "D": "C11H14O2"
    }
    # The provided answer from the LLM is 'B'
    llm_answer_key = "B"
    llm_answer_formula = question_options[llm_answer_key]

    # --- Step 1: Deduce the expected formula from the spectroscopic fragments ---
    # Di-substituted 6-membered aromatic ring -> C6H4
    aromatic_ring = {'C': 6, 'H': 4, 'O': 0}
    
    # Vinyl-H (doublet & doublet of quartets) + one CH3 signal -> propenyl group (-CH=CH-CH3) -> C3H5
    propenyl_group = {'C': 3, 'H': 5, 'O': 0}
    
    # Ester group + second CH3 signal + NO -CH2- signals -> methyl ester (-COOCH3) or acetate (-OCOCH3) -> C2H3O2
    ester_methyl_part = {'C': 2, 'H': 3, 'O': 2}

    # Sum the fragments to get the deduced formula
    deduced_atoms = {}
    for fragment in [aromatic_ring, propenyl_group, ester_methyl_part]:
        for element, count in fragment.items():
            deduced_atoms[element] = deduced_atoms.get(element, 0) + count
            
    # --- Step 2: Compare the LLM's answer with the deduced formula ---
    llm_answer_atoms = parse_formula(llm_answer_formula)
    
    if llm_answer_atoms != deduced_atoms:
        return (f"Incorrect. The formula deduced from the spectral fragments is C{deduced_atoms['C']}H{deduced_atoms['H']}O{deduced_atoms['O']}. "
                f"The provided answer corresponds to {llm_answer_formula}, which does not match.")

    # --- Step 3: Verify constraints to show why the answer is correct and others are wrong ---
    # Constraint 1: Degree of Unsaturation (DoU)
    # Expected DoU = 4 (aromatic ring) + 1 (ester C=O) + 1 (vinyl C=C) = 6
    expected_dou = 6
    
    llm_answer_dou = calculate_dou(llm_answer_atoms)
    if llm_answer_dou != expected_dou:
        return (f"Incorrect. The DoU for the proposed structure should be {expected_dou}. "
                f"The answer {llm_answer_formula} has a DoU of {llm_answer_dou}.")

    # Constraint 2: The "no -CH2- signals" rule, which helps eliminate other options.
    # Let's check option C: C12H14O2, which also has a DoU of 6.
    option_c_atoms = parse_formula(question_options["C"])
    
    # The difference between C12H14O2 and the correct C11H12O2 is C1H2, i.e., a CH2 group.
    diff_c = option_c_atoms.get('C', 0) - deduced_atoms.get('C', 0)
    diff_h = option_c_atoms.get('H', 0) - deduced_atoms.get('H', 0)
    
    if diff_c == 1 and diff_h == 2:
        # This confirms that C12H14O2 contains an extra CH2 unit compared to the structure
        # built from the fragments. This would produce a -CH2- signal, violating the constraint.
        pass # This is expected and confirms our reasoning.
    else:
        # This part of the check is more for validating our own logic.
        # If the difference wasn't CH2, our reasoning for eliminating C12H14O2 would be flawed.
        pass

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check
result = check_answer()
print(result)