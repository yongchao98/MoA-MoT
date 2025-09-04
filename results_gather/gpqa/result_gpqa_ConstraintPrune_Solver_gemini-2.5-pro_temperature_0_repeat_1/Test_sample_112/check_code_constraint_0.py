import re

def check_correctness():
    """
    Checks the correctness of the LLM's answer by verifying it against all
    constraints derived from the spectroscopic data.
    """
    llm_answer_key = "C"
    candidates = {
        "A": "C6H12O2",
        "B": "C6H10O",
        "C": "C6H10O2",
        "D": "C6H12O"
    }
    
    llm_formula = candidates.get(llm_answer_key)
    if not llm_formula:
        return f"The provided answer key '{llm_answer_key}' is not a valid option."

    # Helper function to parse formula and calculate properties
    def get_properties(formula):
        """Parses a formula and returns its properties."""
        props = {'C': 0, 'H': 0, 'O': 0}
        parts = re.findall(r'([A-Z])(\d*)', formula)
        for element, count_str in parts:
            props[element] = int(count_str) if count_str else 1
        
        # Calculate Degree of Unsaturation (DoU)
        dou = props['C'] + 1 - (props['H'] / 2)
        props['DoU'] = dou
        return props

    # --- Constraint 1: Degree of Unsaturation must be >= 2 ---
    # Evidence: C=C (vinyl H) and C=O (1700 cm-1)
    required_dou = 2
    
    # --- Constraint 2: Must be a carboxylic acid ---
    # Evidence: Broad 3000 cm-1 + 1700 cm-1 FTIR peaks, m/z=45 fragment
    # This implies the molecule must have exactly 2 oxygen atoms.
    required_oxygens = 2

    # Check all candidates to find which ones satisfy the constraints
    valid_options = []
    for key, formula in candidates.items():
        props = get_properties(formula)
        
        # Check DoU
        dou_satisfied = props['DoU'] >= required_dou
        
        # Check Oxygen count
        oxygen_satisfied = props['O'] == required_oxygens
        
        if dou_satisfied and oxygen_satisfied:
            valid_options.append(key)

    # Final validation of the LLM's answer
    if llm_answer_key not in valid_options:
        # If the answer is not in the list of valid options, find the reason why it failed.
        props = get_properties(llm_formula)
        if props['DoU'] < required_dou:
            return (f"The answer {llm_formula} ({llm_answer_key}) is incorrect. "
                    f"Its Degree of Unsaturation is {props['DoU']}, but the spectroscopic data "
                    f"requires a DoU of at least {required_dou}.")
        if props['O'] != required_oxygens:
            return (f"The answer {llm_formula} ({llm_answer_key}) is incorrect. "
                    f"It has {props['O']} oxygen atom(s), but the evidence for a carboxylic acid "
                    f"requires {required_oxygens} oxygen atoms.")
        # This case should not be reached if logic is sound
        return f"The answer {llm_formula} ({llm_answer_key}) is incorrect for an unknown reason."

    if len(valid_options) > 1:
        return (f"The answer {llm_answer_key} is plausible but not uniquely correct. "
                f"Options {valid_options} also satisfy all constraints.")

    # If the LLM's answer is the only valid option
    return "Correct"

# Run the check and print the result
result = check_correctness()
print(result)