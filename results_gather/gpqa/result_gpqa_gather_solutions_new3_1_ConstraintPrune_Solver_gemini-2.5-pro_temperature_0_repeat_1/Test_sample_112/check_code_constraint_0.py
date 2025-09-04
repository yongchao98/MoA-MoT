import re

def check_correctness():
    """
    Checks the correctness of the LLM's answer by verifying the molecular formula against
    constraints derived from the spectroscopic data.
    """

    # --- Step 1: Define constraints from the problem description ---
    # Constraint 1: Presence of a carboxylic acid (-COOH) means there must be 2 oxygen atoms.
    required_oxygens = 2
    # Constraint 2: Presence of a C=O and a C=C bond means the Degree of Unsaturation (DoU) must be 2.
    required_dou = 2

    # --- Step 2: Define helper functions ---
    def calculate_dou(c, h):
        """Calculates the Degree of Unsaturation for a formula with C and H."""
        # Formula: DoU = C - H/2 + 1
        return c - (h / 2) + 1

    def parse_formula(formula_str):
        """Parses a chemical formula string (e.g., 'C6H10O2') into a dictionary of atom counts."""
        atoms = {'C': 0, 'H': 0, 'O': 0}
        pattern = r'([A-Z])(\d*)'
        matches = re.findall(pattern, formula_str)
        for element, count in matches:
            # If count is not specified, it's 1
            count = int(count) if count else 1
            if element in atoms:
                atoms[element] += count
        return atoms

    # --- Step 3: Define the problem's options and the answer to check ---
    options = {
        "A": "C6H12O2",
        "B": "C6H10O",
        "C": "C6H10O2",
        "D": "C6H12O"
    }
    
    # The final answer from the LLM to be checked
    llm_answer_text = "<<<C>>>"

    # --- Step 4: Systematically find the correct option ---
    correct_option_letter = None
    for letter, formula in options.items():
        atoms = parse_formula(formula)
        c, h, o = atoms['C'], atoms['H'], atoms['O']

        # Check against Constraint 1 (Oxygen Count)
        if o != required_oxygens:
            continue

        # Check against Constraint 2 (Degree of Unsaturation)
        dou = calculate_dou(c, h)
        if dou != required_dou:
            continue

        # If both constraints are met, this is the correct option
        correct_option_letter = letter
        break

    # --- Step 5: Compare the LLM's answer with the determined correct answer ---
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return f"The provided answer '{llm_answer_text}' is not in the expected format '<<<X>>>'."

    llm_answer_letter = match.group(1)

    if llm_answer_letter == correct_option_letter:
        return "Correct"
    else:
        llm_chosen_formula = options.get(llm_answer_letter, "an invalid option")
        correct_formula = options.get(correct_option_letter, "an unknown correct formula")
        
        # Generate a specific reason why the LLM's choice was wrong
        atoms = parse_formula(llm_chosen_formula)
        o = atoms['O']
        dou = calculate_dou(atoms['C'], atoms['H'])
        
        reason = ""
        if o != required_oxygens:
            reason = f"it has {o} oxygen atoms, but the evidence for a carboxylic acid requires {required_oxygens}."
        elif dou != required_dou:
            reason = f"its Degree of Unsaturation is {int(dou)}, but the evidence for a C=O and a C=C bond requires a DoU of {required_dou}."
        
        return (f"Incorrect. The provided answer is '{llm_answer_letter}', which corresponds to the formula {llm_chosen_formula}. "
                f"The correct answer is '{correct_option_letter}' ({correct_formula}). "
                f"The chosen formula is wrong because {reason}")

# Execute the check and print the result
result = check_correctness()
print(result)