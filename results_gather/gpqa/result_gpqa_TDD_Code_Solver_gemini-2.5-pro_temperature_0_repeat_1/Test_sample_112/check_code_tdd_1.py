import re

def check_answer():
    """
    Checks the correctness of the identified chemical formula based on spectral data.
    """
    llm_answer_option = "C"
    options = {
        "A": "C6H12O",
        "B": "C6H12O2",
        "C": "C6H10O2",
        "D": "C6H10O"
    }

    # --- Step 1: Interpret the spectral data to define constraints ---

    # Constraint 1: Number of Oxygen atoms
    # FTIR broad peak at 3000 cm-1 (O-H stretch) and MS fragment at m/z=45 ([COOH]+)
    # strongly indicate a carboxylic acid group.
    # A carboxylic acid (-COOH) contains 2 oxygen atoms.
    required_oxygens = 2

    # Constraint 2: Degree of Unsaturation (DoU)
    # FTIR peak at 1700 cm-1 indicates a carbonyl group (C=O), which contributes 1 DoU.
    # FTIR peak at 1650 cm-1 and NMR vinyl-hydrogens indicate an alkene (C=C), which contributes 1 DoU.
    # Total required DoU = 1 (for C=O) + 1 (for C=C) = 2.
    required_dou = 2

    # --- Step 2: Check if the LLM's chosen formula satisfies the constraints ---
    
    chosen_formula = options.get(llm_answer_option)
    if not chosen_formula:
        return f"Invalid option '{llm_answer_option}' provided. The option must be one of {list(options.keys())}."

    # Parse the chosen formula
    match = re.match(r'C(\d+)H(\d+)O(\d*)', chosen_formula)
    c = int(match.group(1))
    h = int(match.group(2))
    # If O is not specified, it's 1.
    o = int(match.group(3)) if match.group(3) else 1

    # Check oxygen constraint
    if o != required_oxygens:
        return (f"Incorrect. The answer {chosen_formula} is wrong. "
                f"The spectral data (broad FTIR at 3000 cm-1 and MS peak at m/z=45) indicate a carboxylic acid, "
                f"which requires {required_oxygens} oxygen atoms, but the formula {chosen_formula} has {o}.")

    # Calculate Degree of Unsaturation for the chosen formula
    dou = c - (h / 2) + 1
    
    # Check DoU constraint
    if dou != required_dou:
        return (f"Incorrect. The answer {chosen_formula} is wrong. "
                f"The spectral data (FTIR at 1700 cm-1 and 1650 cm-1, NMR vinyl-H) indicate the presence of both a carbonyl (C=O) and an alkene (C=C), "
                f"requiring a Degree of Unsaturation (DoU) of {required_dou}. "
                f"The formula {chosen_formula} has a DoU of {dou}.")

    # --- Step 3: Verify that no other option also satisfies the constraints ---
    
    valid_options = []
    for option, formula in options.items():
        m = re.match(r'C(\d+)H(\d+)O(\d*)', formula)
        c_i = int(m.group(1))
        h_i = int(m.group(2))
        o_i = int(m.group(3)) if m.group(3) else 1
        dou_i = c_i - (h_i / 2) + 1

        if o_i == required_oxygens and dou_i == required_dou:
            valid_options.append(option)
    
    if len(valid_options) > 1:
        return (f"Ambiguous. The chosen answer {llm_answer_option} is correct, but other options also fit the criteria. "
                f"Options {valid_options} all satisfy the constraints.")

    if llm_answer_option not in valid_options:
        # This case should be caught by the checks above, but it's good practice to include it.
        return f"Incorrect. The correct option should be {valid_options[0]}, but the answer given was {llm_answer_option}."

    return "Correct"

# Run the check
result = check_answer()
print(result)