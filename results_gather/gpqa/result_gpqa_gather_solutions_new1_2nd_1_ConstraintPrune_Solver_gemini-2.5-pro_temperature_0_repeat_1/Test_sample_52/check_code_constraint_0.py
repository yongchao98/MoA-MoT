import re

def parse_formula(formula_str):
    """Parses a chemical formula string into a dictionary of atom counts."""
    pattern = r'([A-Z][a-z]*)(\d*)'
    atoms = {}
    for element, count in re.findall(pattern, formula_str):
        atoms[element] = int(count) if count else 1
    return atoms

def calculate_dou(atoms):
    """Calculates the Degree of Unsaturation for a formula containing C, H, O."""
    c = atoms.get('C', 0)
    h = atoms.get('H', 0)
    # Halogens (X) and Nitrogen (N) are not present, so the formula is simplified.
    # DoU = C + 1 - (H/2) - (X/2) + (N/2)
    return c + 1 - (h / 2)

def check_answer():
    """
    Checks the correctness of the answer based on the problem's constraints.
    """
    # --- Step 1: Deduce the molecular fragments from the spectroscopic data ---

    # "di-substituted 6-membered aromatic ring" -> C6H4
    aromatic_ring = {'C': 6, 'H': 4}
    
    # "FTIR spectrum ... ester group" -> -COO-
    # This accounts for 2 Oxygen atoms and adds 1 to the Degree of Unsaturation (DoU) for the C=O bond.
    ester_group = {'O': 2}
    
    # "1H NMR ... two signals corresponding to vinyl-H (one doublet and one doublet of quartets)"
    # This is a classic signature for a propenyl group: -CH=CH-CH3
    propenyl_group = {'C': 3, 'H': 5} # Adds 1 to DoU for the C=C bond.
    
    # "1H NMR ... two signals corresponding to –CH3 groups"
    # One CH3 is in the propenyl group. The other must be part of the ester.
    
    # "1H NMR ... no signals corresponding to –CH2 groups"
    # This is a critical constraint. To have a second CH3 in the ester without a CH2,
    # it must be a methyl ester: -COOCH3.
    # This fragment adds 1 C (from C=O), 1 C (from O-CH3), and 3 H (from O-CH3).
    # The oxygens are already accounted for.
    methyl_ester_contribution = {'C': 2, 'H': 3}

    # --- Step 2: Assemble the fragments to get the expected formula ---
    
    expected_atoms = {}
    fragments = [aromatic_ring, propenyl_group, methyl_ester_contribution, ester_group]
    for fragment in fragments:
        for atom, count in fragment.items():
            expected_atoms[atom] = expected_atoms.get(atom, 0) + count
            
    expected_formula_str = f"C{expected_atoms['C']}H{expected_atoms['H']}O{expected_atoms['O']}"

    # --- Step 3: Define the options and the given answer ---
    
    options = {
        "A": "C11H14O2",
        "B": "C12H14O2",
        "C": "C11H12O2",
        "D": "C12H12O2"
    }
    
    # The final answer provided by the LLM being checked.
    llm_answer_choice = "C"
    llm_answer_formula = options.get(llm_answer_choice)

    # --- Step 4: Verify the LLM's answer ---

    if not llm_answer_formula:
        return f"Invalid answer choice '{llm_answer_choice}'. The options are A, B, C, D."

    if llm_answer_formula != expected_formula_str:
        return (f"Incorrect. The reasoning leads to the formula {expected_formula_str}, "
                f"but the selected answer corresponds to {llm_answer_formula}.")

    # --- Step 5: Verify that other options are incorrect based on constraints ---
    
    # The required DoU from the structure is: 4 (ring) + 1 (C=C) + 1 (C=O) = 6
    required_dou = 6
    
    for choice, formula_str in options.items():
        if choice == llm_answer_choice:
            continue # Skip the correct answer
            
        atoms = parse_formula(formula_str)
        dou = calculate_dou(atoms)
        
        if dou != required_dou:
            # This check eliminates C11H14O2 (DoU=5) and C12H12O2 (DoU=7)
            continue 
        
        # This check is for C12H14O2, which has the correct DoU of 6.
        # We compare it to the derived correct formula C11H12O2.
        # C12H14O2 = C11H12O2 + CH2
        if atoms['C'] > expected_atoms['C'] and atoms['H'] > expected_atoms['H']:
            if (atoms['C'] - expected_atoms['C'] == 1) and (atoms['H'] - expected_atoms['H'] == 2):
                # This formula contains an extra CH2 group, which is explicitly ruled out.
                continue

        # If an incorrect option passed all checks, our logic is flawed.
        return f"Logic Error: The incorrect option {choice} ({formula_str}) was not properly invalidated."

    return "Correct"

# Run the check
result = check_answer()
print(result)