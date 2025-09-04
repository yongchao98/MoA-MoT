import re

def check_answer_correctness():
    """
    Checks the correctness of the LLM's answer by verifying the given chemical formulas
    against constraints derived from the spectroscopic data.
    """
    # --- Step 1: Define constraints from the problem description ---

    # A broad 3000 cm-1 peak, a 1700 cm-1 peak, and an m/z=45 fragment all point
    # to a carboxylic acid (-COOH), which has 2 oxygen atoms.
    required_oxygens = 2

    # A carboxylic acid contains one C=O bond (1 DoU).
    # A 1650 cm-1 peak and vinyl hydrogens indicate a C=C bond (1 DoU).
    # Total required Degree of Unsaturation (DoU) = 1 + 1 = 2.
    required_dou = 2

    # --- Step 2: Define the options and the LLM's proposed answer ---
    options = {
        "A": "C6H10O",
        "B": "C6H10O2",
        "C": "C6H12O2",
        "D": "C6H12O"
    }
    llm_answer_key = "B"

    # --- Step 3: Define helper functions ---
    def parse_formula(formula_str):
        """Parses a chemical formula string into C, H, O counts."""
        c = int(re.search(r'C(\d+)', formula_str).group(1))
        h = int(re.search(r'H(\d+)', formula_str).group(1))
        o_match = re.search(r'O(\d*)', formula_str)
        if not o_match:
            o = 0
        else:
            # Handles both 'O' (implies O1) and 'O2'
            o = int(o_match.group(1)) if o_match.group(1) else 1
        return c, h, o

    def calculate_dou(c, h):
        """Calculates the Degree of Unsaturation."""
        return c + 1 - (h / 2)

    # --- Step 4: Find the correct option based on constraints ---
    valid_options = []
    for key, formula in options.items():
        try:
            c, h, o = parse_formula(formula)
            dou = calculate_dou(c, h)

            # Check if the formula satisfies both constraints
            if o == required_oxygens and dou == required_dou:
                valid_options.append(key)
        except (AttributeError, ValueError):
            # Handle cases where formula parsing might fail
            continue
    
    # --- Step 5: Validate the LLM's answer ---
    if not valid_options:
        return "Incorrect. None of the provided options satisfy the spectroscopic constraints (2 Oxygens and a DoU of 2)."

    if len(valid_options) > 1:
        return f"Incorrect. The problem is ambiguous as options {valid_options} both satisfy the constraints. The LLM's answer {llm_answer_key} is one of them, but a unique answer cannot be determined."

    correct_key = valid_options[0]
    
    if llm_answer_key == correct_key:
        return "Correct"
    else:
        correct_formula = options[correct_key]
        c, h, o = parse_formula(correct_formula)
        dou = calculate_dou(c, h)
        
        llm_formula = options[llm_answer_key]
        llm_c, llm_h, llm_o = parse_formula(llm_formula)
        llm_dou = calculate_dou(llm_c, llm_h)

        reason = f"Incorrect. The correct answer is {correct_key} ({correct_formula}).\n"
        reason += f"Reasoning: The spectroscopic data requires {required_oxygens} oxygen atoms and a Degree of Unsaturation (DoU) of {required_dou}.\n"
        reason += f"- The LLM's answer {llm_answer_key} ({llm_formula}) has {llm_o} oxygen(s) and a DoU of {llm_dou}.\n"
        reason += f"- The correct answer {correct_key} ({correct_formula}) has {o} oxygen(s) and a DoU of {dou}, which matches all constraints."
        return reason

# Execute the check and print the result
result = check_answer_correctness()
print(result)