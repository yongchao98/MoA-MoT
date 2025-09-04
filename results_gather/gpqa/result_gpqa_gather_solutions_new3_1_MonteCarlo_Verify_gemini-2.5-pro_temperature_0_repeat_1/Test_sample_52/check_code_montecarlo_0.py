import re

def check_answer():
    """
    Checks the correctness of the final answer for the given chemistry problem.
    """
    
    # --- Define Constraints from the Question ---

    # 1. Degree of Unsaturation (DoU) from structural features
    # Aromatic ring = 4 (1 ring + 3 double bonds)
    # Ester group (C=O) = 1
    # Vinyl group (C=C) = 1
    required_dou = 4 + 1 + 1

    # 2. Atom count from NMR fragments
    # Aromatic ring: di-substituted -> C6H4
    # Vinyl group: -CH=CH-CH3 (propenyl) -> C3H5
    # Ester + second methyl group (with no -CH2-): -COOCH3 or -OCOCH3 -> C2H3O2
    expected_atoms = {
        'C': 6 + 3 + 2,
        'H': 4 + 5 + 3,
        'O': 2
    }
    
    # --- The Answer to be Checked ---
    
    # The question options are:
    # A) C11H14O2, B) C12H14O2, C) C11H12O2, D) C12H12O2
    # The final answer provided is <<<C>>>, which corresponds to C11H12O2.
    final_answer_formula = "C11H12O2"

    # --- Helper Functions ---

    def parse_formula(formula_str):
        """Parses a chemical formula string into a dictionary of atom counts."""
        atoms = {'C': 0, 'H': 0, 'O': 0}
        # Find all element-count pairs (e.g., 'C11', 'H12', 'O2')
        parts = re.findall(r'([A-Z][a-z]*)(\d*)', formula_str)
        for element, count in parts:
            if element in atoms:
                atoms[element] = int(count) if count else 1
        return atoms

    def calculate_dou(atoms):
        """Calculates the Degree of Unsaturation for a given set of atoms."""
        # Formula: DoU = C + 1 - (H/2) - (X/2) + (N/2)
        # For C, H, O: DoU = C + 1 - H/2
        if 'C' not in atoms or 'H' not in atoms:
            return None
        return atoms.get('C', 0) + 1 - (atoms.get('H', 0) / 2)

    # --- Verification Logic ---

    try:
        # 1. Parse the formula of the chosen answer
        answer_atoms = parse_formula(final_answer_formula)

        # 2. Check Oxygen count (basic sanity check)
        if answer_atoms.get('O', 0) != 2:
            return f"Incorrect: The answer {final_answer_formula} has {answer_atoms.get('O', 0)} oxygen atoms, but the presence of an ester group requires 2 oxygen atoms."

        # 3. Check Degree of Unsaturation
        dou = calculate_dou(answer_atoms)
        if dou != required_dou:
            return f"Incorrect: The answer {final_answer_formula} has a Degree of Unsaturation (DoU) of {dou}. However, the structural features (aromatic ring, ester, vinyl group) require a DoU of {required_dou}."

        # 4. Check atom count based on NMR fragment analysis
        # This is the most specific check
        if answer_atoms['C'] != expected_atoms['C']:
            return f"Incorrect: The answer {final_answer_formula} has {answer_atoms['C']} carbon atoms. However, summing the fragments from the NMR data (C6H4 ring + C3H5 propenyl group + C2H3O2 ester/methyl part) results in {expected_atoms['C']} carbon atoms."
        
        if answer_atoms['H'] != expected_atoms['H']:
            return f"Incorrect: The answer {final_answer_formula} has {answer_atoms['H']} hydrogen atoms. However, summing the fragments from the NMR data (C6H4 ring + C3H5 propenyl group + C2H3O2 ester/methyl part) results in {expected_atoms['H']} hydrogen atoms."

        # If all checks pass, the answer is correct.
        # We can also add a check for why other options are wrong.
        # For example, C12H14O2 has the correct DoU but differs by a CH2 group, which is explicitly forbidden.
        # The chosen answer C11H12O2 is the only one that satisfies all constraints.
        
        return "Correct"

    except Exception as e:
        return f"An error occurred during verification: {e}"

# Run the check and print the result
result = check_answer()
print(result)