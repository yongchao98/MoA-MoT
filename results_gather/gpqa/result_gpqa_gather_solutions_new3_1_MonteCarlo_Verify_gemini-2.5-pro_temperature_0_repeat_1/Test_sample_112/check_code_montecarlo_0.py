import re

def check_correctness():
    """
    Checks the correctness of the answer to the chemical identification problem.
    """
    # 1. Define constraints based on the spectroscopic data from the question.
    # - Carboxylic acid (-COOH) is indicated by:
    #   - FTIR: very broad peak at 3000 cm⁻¹, strong peak at 1700 cm⁻¹
    #   - Mass Spec: fragment at m/z = 45 ([COOH]⁺)
    # - Alkene (C=C) is indicated by:
    #   - FTIR: strong peak at 1650 cm⁻¹
    #   - 1H NMR: peaks for vinyl-hydrogens
    
    # A carboxylic acid requires 2 oxygen atoms and has 1 degree of unsaturation (DoU) from C=O.
    # An alkene has 1 degree of unsaturation from C=C.
    required_oxygens = 2
    required_dou = 1 + 1  # 1 for C=O, 1 for C=C

    # 2. Define the candidate options and the provided answer.
    options = {
        "A": "C6H10O",
        "B": "C6H12O",
        "C": "C6H12O2",
        "D": "C6H10O2"
    }
    provided_answer_letter = "D"

    # 3. Helper functions to parse formulas and calculate DoU.
    def parse_formula(formula_str):
        """Parses a chemical formula string into a dictionary of atom counts."""
        atoms = {'C': 0, 'H': 0, 'O': 0}
        # This regex finds atom-count pairs (e.g., C6, H12, O)
        pattern = r'([A-Z])(\d*)'
        for element, count in re.findall(pattern, formula_str):
            if element in atoms:
                # If count is empty, it means 1 atom (e.g., 'O' in C6H10O)
                atoms[element] = int(count) if count else 1
        return atoms

    def calculate_dou(c, h):
        """Calculates the Degree of Unsaturation for a formula with C, H, O."""
        return c - h / 2 + 1

    # 4. Evaluate each option against the constraints.
    correct_option = None
    failure_reasons = {}

    for letter, formula in options.items():
        atom_counts = parse_formula(formula)
        c, h, o = atom_counts['C'], atom_counts['H'], atom_counts['O']

        # Check 1: Number of oxygen atoms
        if o != required_oxygens:
            failure_reasons[letter] = f"Formula {formula} has {o} oxygen atom(s), but the evidence for a carboxylic acid requires {required_oxygens}."
            continue
            
        # Check 2: Degree of Unsaturation (DoU)
        dou = calculate_dou(c, h)
        if dou != required_dou:
            failure_reasons[letter] = f"Formula {formula} has a Degree of Unsaturation of {int(dou)}, but the evidence for both a C=O and a C=C bond requires a DoU of {required_dou}."
            continue
            
        # If a formula passes all checks, it's the correct one.
        correct_option = letter

    # 5. Compare the identified correct option with the provided answer.
    if provided_answer_letter == correct_option:
        return "Correct"
    else:
        reason = failure_reasons.get(provided_answer_letter, "The provided answer does not satisfy the problem's constraints.")
        return f"Incorrect. The provided answer '{provided_answer_letter}' is wrong. {reason}"

# Execute the check and print the result.
result = check_correctness()
print(result)