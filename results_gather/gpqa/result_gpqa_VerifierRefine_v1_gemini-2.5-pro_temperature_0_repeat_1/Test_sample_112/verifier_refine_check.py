import re

def check_answer():
    """
    This function checks the correctness of the given answer by verifying it against
    the constraints derived from the spectroscopic data in the question.
    """
    # --- 1. Define Constraints from the Question ---

    # From FTIR:
    # - Broad peak at 3000 cm⁻¹ -> O-H stretch of a carboxylic acid.
    # - Strong peak at 1700 cm⁻¹ -> C=O stretch, consistent with a carboxylic acid.
    # - Strong peak at 1650 cm⁻¹ -> C=C stretch of an alkene.
    # From 1H NMR:
    # - Vinyl-hydrogens -> Confirms the presence of a C=C double bond.
    # From Mass Spectrometry:
    # - Fragment at m/z = 45 -> Corresponds to the [COOH]⁺ fragment.

    # Synthesis of constraints:
    # Constraint A: The compound must be a carboxylic acid. A carboxylic acid group (-COOH)
    # contains exactly two oxygen atoms.
    required_oxygens = 2

    # Constraint B: The compound contains a C=C double bond (alkene) in addition to the
    # C=O double bond of the carboxylic acid. This means the molecule must have at least
    # two degrees of unsaturation (one for C=O, one for C=C).
    required_dou = 2

    # --- 2. Define the Given Answer and Options ---
    given_answer_key = "B"
    options = {
        "A": "C6H10O",
        "B": "C6H10O2",
        "C": "C6H12O",
        "D": "C6H12O2"
    }
    
    answer_formula = options.get(given_answer_key)

    # --- 3. Helper Functions ---
    def parse_formula(formula):
        """Parses a chemical formula string into a dictionary of element counts."""
        elements = {'C': 0, 'H': 0, 'O': 0}
        # Finds all element-count pairs (e.g., C6, H12, O2)
        parts = re.findall(r'([A-Z][a-z]*)(\d*)', formula)
        for element, count in parts:
            if element in elements:
                elements[element] = int(count) if count else 1
        return elements

    def calculate_dou(counts):
        """Calculates the Degree of Unsaturation for a CxHyOz compound."""
        C = counts.get('C', 0)
        H = counts.get('H', 0)
        # Formula: DoU = C + 1 - (H/2)
        return C + 1 - (H / 2)

    # --- 4. Perform Checks ---
    
    # Parse the formula of the given answer
    try:
        counts = parse_formula(answer_formula)
    except Exception as e:
        return f"Error: Could not parse the formula '{answer_formula}': {e}"

    # Check Constraint A: Number of Oxygen Atoms
    num_oxygens = counts.get('O', 0)
    if num_oxygens != required_oxygens:
        return (f"Incorrect. The answer {answer_formula} has {num_oxygens} oxygen atom(s). "
                f"The combined FTIR and Mass Spec data (broad O-H, C=O, and m/z=45 fragment) "
                f"strongly indicate a carboxylic acid, which requires {required_oxygens} oxygen atoms.")

    # Check Constraint B: Degree of Unsaturation (DoU)
    dou = calculate_dou(counts)
    if dou != required_dou:
        return (f"Incorrect. The answer {answer_formula} has a Degree of Unsaturation (DoU) of {int(dou)}. "
                f"The data indicates the presence of both a carbonyl group (C=O) and an alkene (C=C), "
                f"which requires a total DoU of {required_dou}.")

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check and print the result
result = check_answer()
print(result)