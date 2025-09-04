import re

def check_correctness():
    """
    This function checks the correctness of the LLM's answer to the chemistry question.
    It derives constraints from the spectroscopic data and verifies which of the given
    chemical formulas satisfies all constraints.
    """
    
    # 1. Define the problem space from the question
    # The options as listed in the question:
    # A) C6H10O, B) C6H12O, C) C6H12O2, D) C6H10O2
    options = {
        "A": "C6H10O",
        "B": "C6H12O",
        "C": "C6H12O2",
        "D": "C6H10O2"
    }
    
    # The final answer provided by the LLM to be checked
    llm_answer_letter = "D"

    # 2. Derive constraints from the spectroscopic data provided in the question
    # - FTIR: Broad peak at 3000 cm-1, strong peak at 1700 cm-1 -> Carboxylic acid (-COOH)
    # - FTIR: Strong peak at 1650 cm-1 -> Alkene (C=C)
    # - 1H NMR: Vinyl-hydrogens -> Confirms Alkene (C=C)
    # - Mass Spec: Fragment at m/z = 45 -> Confirms Carboxylic acid ([COOH]+)
    
    # Constraint 1: The presence of a carboxylic acid group (-COOH) means the molecule must have 2 oxygen atoms.
    required_oxygens = 2
    
    # Constraint 2: The Degree of Unsaturation (DoU) must account for both the C=O and C=C bonds.
    # DoU from C=O in carboxylic acid = 1
    # DoU from C=C in alkene = 1
    # Total required DoU = 1 + 1 = 2
    required_dou = 2

    # 3. Define helper functions for analysis
    def parse_formula(formula_str):
        """Parses a chemical formula string to get the count of C, H, and O atoms."""
        c_match = re.search(r'C(\d+)', formula_str)
        h_match = re.search(r'H(\d+)', formula_str)
        o_match = re.search(r'O(\d*)', formula_str)
        
        carbons = int(c_match.group(1)) if c_match else 0
        hydrogens = int(h_match.group(1)) if h_match else 0
        
        if o_match:
            # Handles cases like 'O' (1 oxygen) and 'O2' (2 oxygens)
            oxygens_str = o_match.group(1)
            oxygens = int(oxygens_str) if oxygens_str else 1
        else:
            oxygens = 0
            
        return carbons, hydrogens, oxygens

    def calculate_dou(carbons, hydrogens):
        """Calculates the Degree of Unsaturation using the formula: C + 1 - (H/2)."""
        return carbons + 1 - (hydrogens / 2)

    # 4. Evaluate all options against the constraints to find the scientifically correct answer
    scientifically_correct_letter = None
    for letter, formula in options.items():
        c, h, o = parse_formula(formula)
        dou = calculate_dou(c, h)
        
        # Check if both constraints are met
        if o == required_oxygens and dou == required_dou:
            scientifically_correct_letter = letter
            break # Assuming only one correct answer exists

    # 5. Compare the LLM's answer with the derived correct answer and provide a verdict
    if llm_answer_letter == scientifically_correct_letter:
        return "Correct"
    else:
        # Analyze the LLM's chosen formula to explain why it's wrong
        llm_formula = options.get(llm_answer_letter)
        if not llm_formula:
            return f"Incorrect. The provided answer '{llm_answer_letter}' is not one of the valid options."

        c, h, o = parse_formula(llm_formula)
        dou = calculate_dou(c, h)
        
        reasons = []
        if o != required_oxygens:
            reasons.append(f"its formula ({llm_formula}) has {o} oxygen atom(s), but the spectroscopic data indicates a carboxylic acid, which requires {required_oxygens} oxygen atoms")
        if dou != required_dou:
            reasons.append(f"its formula ({llm_formula}) has a Degree of Unsaturation of {int(dou)}, but the data indicates both a C=O and a C=C bond, requiring a Degree of Unsaturation of {required_dou}")
        
        reason_str = " and ".join(reasons)
        correct_formula = options[scientifically_correct_letter]
        
        return f"Incorrect. The provided answer '{llm_answer_letter}' is wrong because {reason_str}. The correct answer is '{scientifically_correct_letter}' ({correct_formula}), as it is the only option with 2 oxygens and a Degree of Unsaturation of 2."

# Execute the check and print the result
result = check_correctness()
print(result)