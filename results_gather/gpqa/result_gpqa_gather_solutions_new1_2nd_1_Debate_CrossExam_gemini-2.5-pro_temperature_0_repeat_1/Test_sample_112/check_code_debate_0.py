import re

def check_correctness():
    """
    This function checks the correctness of the provided LLM answer based on the chemical analysis of the question.
    
    The analysis is as follows:
    1.  **FTIR Data**:
        - Broad peak at 3000 cm⁻¹ + strong peak at 1700 cm⁻¹ -> Carboxylic acid (-COOH).
        - Strong peak at 1650 cm⁻¹ -> Alkene (C=C).
    2.  **1H NMR Data**:
        - Vinyl-hydrogens -> Confirms alkene (C=C).
    3.  **Mass Spectrometry Data**:
        - Fragment at m/z = 45 -> Confirms carboxylic acid ([COOH]⁺).

    **Constraints derived from the data**:
    1.  The molecule must be a carboxylic acid, which means it must contain **2 oxygen atoms**.
    2.  The molecule must contain both a C=O double bond (from the acid) and a C=C double bond (from the alkene). This requires a total **Degree of Unsaturation (DoU) of 2**.
    """
    
    # The final answer provided by the LLM
    llm_answer_text = "<<<C>>>"
    
    # The options as defined in the question
    options = {
        'A': 'C6H12O',
        'B': 'C6H10O',
        'C': 'C6H10O2',
        'D': 'C6H12O2'
    }
    
    # Extract the letter choice from the LLM's answer
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Incorrect: The answer format is invalid. It should be in the format '<<<X>>>' where X is one of the options A, B, C, or D."
    
    chosen_option = match.group(1)
    chosen_formula = options[chosen_option]

    # Helper function to parse a chemical formula string
    def parse_formula(formula):
        """Parses a formula string into a dictionary of element counts."""
        counts = {}
        # Regex to find elements and their optional counts
        pattern = r'([A-Z][a-z]*)(\d*)'
        for element, count in re.findall(pattern, formula):
            counts[element] = int(count) if count else 1
        return counts

    # Helper function to calculate the Degree of Unsaturation (DoU)
    def calculate_dou(counts):
        """Calculates DoU for a formula containing C, H, O."""
        # Formula: DoU = C + 1 - (H/2)
        # Halogens and Nitrogen are not present, so the formula is simpler.
        c = counts.get('C', 0)
        h = counts.get('H', 0)
        return c + 1 - (h / 2)

    # --- Verification Steps ---

    # 1. Verify the number of oxygen atoms
    formula_counts = parse_formula(chosen_formula)
    oxygen_count = formula_counts.get('O', 0)
    
    if oxygen_count != 2:
        return (f"Incorrect: The chosen formula {chosen_formula} (Option {chosen_option}) has {oxygen_count} oxygen atom(s). "
                "The evidence from FTIR (broad peak at 3000 cm⁻¹) and Mass Spectrometry (fragment at m/z=45) "
                "strongly indicates a carboxylic acid, which must contain 2 oxygen atoms.")

    # 2. Verify the Degree of Unsaturation (DoU)
    dou = calculate_dou(formula_counts)
    
    if dou != 2:
        return (f"Incorrect: The chosen formula {chosen_formula} (Option {chosen_option}) has a Degree of Unsaturation (DoU) of {int(dou)}. "
                "The evidence indicates the presence of both a carboxylic acid (C=O bond, 1 DoU) and an alkene (C=C bond, 1 DoU), "
                "requiring a total DoU of 2.")

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check and print the result
result = check_correctness()
print(result)