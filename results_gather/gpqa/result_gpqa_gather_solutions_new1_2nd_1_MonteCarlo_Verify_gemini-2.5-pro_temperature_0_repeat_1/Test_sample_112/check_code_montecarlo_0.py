import re

def check_answer():
    """
    Checks the correctness of the LLM's answer based on the chemical constraints.
    """
    # --- Problem Setup ---
    # The question provides spectroscopic data and four possible chemical formulas.
    # The final answer block maps the formulas to letters A, B, C, D.
    # We will use the mapping from the final provided answer block for consistency.
    question_options = {
        "A": "C6H10O2",
        "B": "C6H10O",
        "C": "C6H12O2",
        "D": "C6H12O"
    }
    llm_answer_letter = "A"

    # --- Deriving Constraints from Spectroscopic Data ---
    # 1. FTIR (broad 3000, 1700) + Mass Spec (m/z=45) => Carboxylic Acid (-COOH)
    #    - This means the molecule must have exactly 2 oxygen atoms.
    # 2. FTIR (1650) + 1H NMR (vinyl-H) => Alkene (C=C)
    # 3. The carboxylic acid has a C=O bond.
    # 4. Degree of Unsaturation (DoU) = (number of C=O bonds) + (number of C=C bonds)
    #    - DoU must be at least 1 (for C=O) + 1 (for C=C) = 2.
    #    - Assuming no rings, the DoU must be exactly 2.
    
    constraints = {
        "oxygen_count": 2,
        "degree_of_unsaturation": 2
    }

    # --- Helper Functions ---
    def parse_formula(formula):
        """Parses a chemical formula string into a dictionary of atom counts."""
        atoms = {'C': 0, 'H': 0, 'O': 0}
        # Find all atom-count pairs (e.g., C6, H10, O2)
        parts = re.findall(r'([A-Z][a-z]*)(\d*)', formula)
        for atom, count in parts:
            if atom in atoms:
                # If count is empty, it's 1. Otherwise, convert to int.
                atoms[atom] = int(count) if count else 1
        return atoms

    def calculate_dou(atoms):
        """Calculates the Degree of Unsaturation for a formula CxHyOz."""
        # Formula: DoU = C + 1 - (H/2) - (X/2) + (N/2)
        # For CHO compounds, it simplifies to:
        C = atoms.get('C', 0)
        H = atoms.get('H', 0)
        return C + 1 - (H / 2)

    # --- Verification Logic ---
    correct_options = []
    error_messages = []

    for option_letter, formula in question_options.items():
        atoms = parse_formula(formula)
        oxygen_count = atoms.get('O', 0)
        dou = calculate_dou(atoms)
        
        is_correct = True
        reasons = []
        
        # Check Oxygen Count
        if oxygen_count != constraints["oxygen_count"]:
            is_correct = False
            reasons.append(f"has {oxygen_count} oxygen(s), but the evidence for a carboxylic acid requires {constraints['oxygen_count']}.")
        
        # Check Degree of Unsaturation
        if dou != constraints["degree_of_unsaturation"]:
            is_correct = False
            reasons.append(f"has a Degree of Unsaturation (DoU) of {dou}, but the evidence for both a C=O and a C=C bond requires a DoU of {constraints['degree_of_unsaturation']}.")

        if is_correct:
            correct_options.append(option_letter)
        else:
            error_messages.append(f"Option {option_letter} ({formula}) is incorrect because it " + " and ".join(reasons))

    # --- Final Decision ---
    if len(correct_options) != 1:
        return f"Analysis Error: Found {len(correct_options)} options that satisfy all constraints. Correct options: {correct_options}. Full analysis:\n" + "\n".join(error_messages)

    actual_correct_letter = correct_options[0]

    if llm_answer_letter == actual_correct_letter:
        return "Correct"
    else:
        return (f"Incorrect. The provided answer is {llm_answer_letter}, but the only formula that satisfies all constraints is {question_options[actual_correct_letter]} (Option {actual_correct_letter}).\n\n"
                f"Reasoning:\n"
                f"1. The evidence (FTIR, Mass Spec) points to a carboxylic acid, which requires 2 oxygen atoms.\n"
                f"2. The evidence (FTIR, NMR) points to an alkene (C=C bond) in addition to the acid's C=O bond, requiring a total Degree of Unsaturation of 2.\n\n"
                f"Analysis of the provided answer ({llm_answer_letter}):\n"
                f"{[msg for msg in error_messages if msg.startswith(f'Option {llm_answer_letter}')][0]}\n\n"
                f"Analysis of all options:\n" + "\n".join(error_messages))

# Run the check
result = check_answer()
print(result)