import re

def check_correctness():
    """
    This function checks the correctness of the provided answer based on the constraints
    derived from the spectroscopic data in the question.
    """
    
    # --- Step 1: Define constraints from the problem description ---
    # Evidence for carboxylic acid (-COOH):
    # - FTIR: Very broad peak at 3000 cm⁻¹ (O-H) and strong peak at 1700 cm⁻¹ (C=O).
    # - Mass Spec: Fragment at m/z = 45 ([COOH]⁺).
    # This implies the formula must have 2 oxygen atoms.
    required_oxygens = 2
    
    # Evidence for degrees of unsaturation (DoU):
    # - One C=O double bond from the carboxylic acid (1 DoU).
    # - One C=C double bond from the alkene group (FTIR at 1650 cm⁻¹, NMR vinyl-H) (1 DoU).
    # Total required DoU is 2.
    required_dou = 2

    # --- Step 2: Define the options and the provided answer ---
    options = {
        'A': 'C6H10O',
        'B': 'C6H10O2',
        'C': 'C6H12O',
        'D': 'C6H12O2'
    }
    # The final answer provided by the LLM.
    llm_answer = 'B'

    # --- Step 3: Define helper functions ---
    def parse_formula(formula_str):
        """Parses a chemical formula string into a dictionary of atom counts."""
        counts = {'C': 0, 'H': 0, 'O': 0}
        # Use regex to find all element-count pairs (e.g., C6, H10, O2)
        parts = re.findall(r'([A-Z])(\d*)', formula_str)
        for element, count_str in parts:
            count = int(count_str) if count_str else 1
            if element in counts:
                counts[element] = count
        return counts

    def calculate_dou(counts):
        """Calculates the Degree of Unsaturation for a C, H, O formula."""
        # Formula: DoU = C + 1 - (H/2)
        C = counts.get('C', 0)
        H = counts.get('H', 0)
        return C + 1 - (H / 2)

    # --- Step 4: Analyze each option to find the correct one ---
    correct_option = None
    for option_letter, formula in options.items():
        counts = parse_formula(formula)
        
        # Check if the number of oxygens matches the requirement
        oxygen_check = (counts.get('O', 0) == required_oxygens)
        
        # Check if the Degree of Unsaturation matches the requirement
        dou_check = (calculate_dou(counts) == required_dou)
        
        if oxygen_check and dou_check:
            correct_option = option_letter
            break # Found the unique correct option

    # --- Step 5: Compare the LLM's answer with the derived correct answer ---
    if correct_option is None:
        return "Analysis Error: No option satisfies the constraints derived from the spectral data."

    if llm_answer == correct_option:
        return "Correct"
    else:
        # Provide a reason why the LLM's answer is wrong.
        incorrect_formula = options[llm_answer]
        incorrect_counts = parse_formula(incorrect_formula)
        incorrect_oxygens = incorrect_counts.get('O', 0)
        incorrect_dou = calculate_dou(incorrect_counts)
        
        reason = f"Incorrect. The provided answer is {llm_answer} ({incorrect_formula}). "
        
        if incorrect_oxygens != required_oxygens:
            reason += (f"This formula has {incorrect_oxygens} oxygen atom(s), but the evidence for a carboxylic acid "
                       f"requires {required_oxygens} oxygen atoms. ")
        
        if incorrect_dou != required_dou:
            reason += (f"This formula has a Degree of Unsaturation (DoU) of {incorrect_dou}, but the evidence for "
                       f"both a C=O and a C=C double bond requires a DoU of {required_dou}. ")
            
        reason += f"The correct option is {correct_option} ({options[correct_option]}), which satisfies all constraints."
        return reason

# Execute the check and print the result
print(check_correctness())