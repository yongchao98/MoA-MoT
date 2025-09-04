import re

def check_correctness():
    """
    Checks the correctness of the LLM's answer based on the provided chemical analysis question.
    """
    # The final answer provided by the LLM.
    llm_answer_text = "<<<C>>>"

    # Define the options as provided in the question.
    options = {
        "A": "C6H10O",
        "B": "C6H12O",
        "C": "C6H10O2",
        "D": "C6H12O2"
    }

    # --- Step 1: Parse the LLM's answer ---
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Incorrect: The answer format is invalid. It should be in the format <<<X>>> where X is A, B, C, or D."
    
    llm_choice_letter = match.group(1)
    llm_chosen_formula = options.get(llm_choice_letter)

    if not llm_chosen_formula:
        return f"Incorrect: The chosen option '{llm_choice_letter}' is not a valid option."

    # --- Step 2: Define constraints from the spectroscopic data ---
    # Evidence: FTIR (broad 3000, 1700, 1650 cm-1), NMR (vinyl-H), MS (m/z=45)
    # Conclusion: The molecule is an unsaturated carboxylic acid.
    
    # Constraint 1: A carboxylic acid (-COOH) requires 2 oxygen atoms.
    required_oxygens = 2
    
    # Constraint 2: A carboxylic acid (C=O) and an alkene (C=C) require a Degree of Unsaturation (DoU) of 2.
    required_dou = 2

    # --- Step 3: Helper functions to analyze formulas ---
    def parse_formula(formula_str):
        """Extracts the number of C, H, and O atoms from a formula string."""
        c = int(re.search(r'C(\d+)', formula_str).group(1))
        h = int(re.search(r'H(\d+)', formula_str).group(1))
        o_match = re.search(r'O(\d*)', formula_str)
        o = int(o_match.group(1)) if o_match.group(1) and o_match.group(1).isdigit() else 1
        return c, h, o

    def calculate_dou(c, h):
        """Calculates the Degree of Unsaturation."""
        return c - (h / 2) + 1

    # --- Step 4: Evaluate all candidate formulas to find the correct one ---
    correct_formulas = []
    for formula in options.values():
        c, h, o = parse_formula(formula)
        
        # Check if the formula meets both constraints
        if o == required_oxygens and calculate_dou(c, h) == required_dou:
            correct_formulas.append(formula)

    # There should be only one correct formula among the options
    if len(correct_formulas) != 1:
        return f"Analysis Error: The problem is ambiguous. Found {len(correct_formulas)} formulas that fit the criteria: {correct_formulas}."
        
    correct_formula = correct_formulas[0]

    # --- Step 5: Compare the LLM's answer with the derived correct answer ---
    if llm_chosen_formula == correct_formula:
        return "Correct"
    else:
        # Generate a detailed reason for the incorrectness
        llm_c, llm_h, llm_o = parse_formula(llm_chosen_formula)
        llm_dou = calculate_dou(llm_c, llm_h)
        
        reason = f"The chosen formula {llm_chosen_formula} is incorrect. "
        if llm_o != required_oxygens:
            reason += f"It has {llm_o} oxygen atom(s), but the spectroscopic data (FTIR broad O-H, C=O peak, and MS m/z=45 fragment) indicate a carboxylic acid, which requires {required_oxygens} oxygen atoms. "
        if llm_dou != required_dou:
            reason += f"It has a Degree of Unsaturation (DoU) of {int(llm_dou)}, but the data (C=O from the acid and C=C from the alkene) requires a DoU of {required_dou}. "
        
        reason += f"The correct formula is {correct_formula}."
        return f"Incorrect: {reason}"

# Run the check
result = check_correctness()
print(result)