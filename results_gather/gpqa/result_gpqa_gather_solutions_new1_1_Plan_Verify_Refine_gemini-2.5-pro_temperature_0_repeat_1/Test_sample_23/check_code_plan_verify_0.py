import re

def check_correctness():
    """
    This function checks the correctness of the provided answer by verifying each constraint
    of the chemistry problem.
    """

    # --- Define Chemical Properties (Simulated Knowledge Base) ---
    # This dictionary stores the properties of the relevant compounds.
    properties = {
        'C6H12': {'name': 'cyclohexane', 'saturated': True, 'aromatic': False, 'conjugated': False, 'decolorizes_br2': False},
        'C6H6': {'name': 'benzene', 'saturated': False, 'aromatic': True, 'conjugated': True, 'decolorizes_br2': False},
        'C6H10': {'name': 'cyclohexene', 'saturated': False, 'aromatic': False, 'conjugated': False, 'decolorizes_br2': True},
        'C6H8': {'name': '1,4-cyclohexadiene', 'saturated': False, 'aromatic': False, 'conjugated': False, 'decolorizes_br2': True},
        'C6H8_conjugated': {'name': '1,3-cyclohexadiene', 'saturated': False, 'aromatic': False, 'conjugated': True, 'decolorizes_br2': True}
    }

    def parse_formula(formula):
        """Parses a chemical formula like 'C6H12' into a dictionary {'C': 6, 'H': 12}."""
        parts = re.findall(r'([A-Z][a-z]*)(\d*)', formula)
        atoms = {}
        for element, count in parts:
            atoms[element] = int(count) if count else 1
        return atoms

    def calculate_h_mass_fraction(formula):
        """Calculates the mass fraction of hydrogen in a hydrocarbon using atomic masses C=12, H=1."""
        atoms = parse_formula(formula)
        mass_c = atoms.get('C', 0) * 12
        mass_h = atoms.get('H', 0) * 1
        total_mass = mass_c + mass_h
        return mass_h / total_mass if total_mass > 0 else 0

    # --- Step 1: Verify Substance Z ---
    # Hypothesis from the reasoning: Z is cyclohexane (C6H12).
    Z_formula = 'C6H12'
    
    # Constraint: Mass fraction of H is 14.28% (~1/7).
    h_fraction = calculate_h_mass_fraction(Z_formula)
    if not (0.142 < h_fraction < 0.143):
        return f"Incorrect: The mass fraction of hydrogen for the proposed Z ({Z_formula}) is {h_fraction:.4f}, which is not approximately 14.28% (1/7)."
    
    # Constraint: Z is saturated (does not react further with hydrogen).
    if not properties[Z_formula]['saturated']:
        return f"Incorrect: The proposed Z ({Z_formula}) is not saturated, but the problem states it is."

    # --- Step 2: Verify Mixture Y ---
    # Hypothesis from the reasoning: Y = {cyclohexane (C6H12), benzene (C6H6)}.
    Y_formulas = {'C6H12', 'C6H6'}

    # Constraint: Z is a constituent of Y.
    if Z_formula not in Y_formulas:
        return f"Incorrect: The proposed mixture Y {Y_formulas} does not contain Z ({Z_formula})."
    
    # Constraint: Y does not decolorize bromine water.
    for formula in Y_formulas:
        if properties[formula]['decolorizes_br2']:
            return f"Incorrect: A component of mixture Y ({formula}) decolorizes bromine water, which violates the problem constraints."
    
    # Constraint: Hydrogenation of Y gives only Z.
    # Hydrogenation of C6H12 -> C6H12 (Z)
    # Hydrogenation of C6H6 -> C6H12 (Z)
    # This condition is met by the hypothesis.

    # --- Step 3: Verify Mixture X ---
    # Hypothesis from the reasoning: X = {cyclohexene (C6H10), 1,4-cyclohexadiene (C6H8)}.
    X_formulas = {'C6H10', 'C6H8'}

    # Constraint: X decolorizes bromine water.
    for formula in X_formulas:
        if not properties[formula]['decolorizes_br2']:
            return f"Incorrect: A component of mixture X ({formula}) does not decolorize bromine water, which violates the problem constraints."
    
    # Constraint: No conjugated multiple bonds in X.
    for formula in X_formulas:
        if properties[formula]['conjugated']:
            return f"Incorrect: A component of mixture X ({formula}) has conjugated bonds, which violates the problem constraints."

    # Constraint: Hydrogenation of X gives only Z.
    # Hydrogenation of C6H10 -> C6H12 (Z)
    # Hydrogenation of C6H8 -> C6H12 (Z)
    # This condition is met by the hypothesis.

    # --- Step 4: Verify the Disproportionation Reaction (Atom Conservation) ---
    # Constraint: X -> Y must be a balanced reaction.
    atoms_X = {'C': 0, 'H': 0}
    for formula in X_formulas:
        parsed = parse_formula(formula)
        atoms_X['C'] += parsed.get('C', 0)
        atoms_X['H'] += parsed.get('H', 0)

    atoms_Y = {'C': 0, 'H': 0}
    for formula in Y_formulas:
        parsed = parse_formula(formula)
        atoms_Y['C'] += parsed.get('C', 0)
        atoms_Y['H'] += parsed.get('H', 0)

    if atoms_X != atoms_Y:
        return f"Incorrect: The disproportionation reaction is not balanced. Atoms in X: {atoms_X}, Atoms in Y: {atoms_Y}."

    # --- Step 5: Final Calculation and Verification ---
    # The question asks for the total number of hydrogen atoms in the two liquids of mixture X.
    total_h_atoms_calculated = atoms_X['H']
    
    # The provided answer is 'B', which corresponds to 18.
    final_answer_value = 18

    if total_h_atoms_calculated != final_answer_value:
        return f"Incorrect: The calculated total number of hydrogen atoms in mixture X is {total_h_atoms_calculated}, but the provided answer corresponds to {final_answer_value}."

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check
result = check_correctness()
print(result)