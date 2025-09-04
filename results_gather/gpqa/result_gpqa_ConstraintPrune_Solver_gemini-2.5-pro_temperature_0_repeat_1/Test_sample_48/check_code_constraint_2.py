import re
from collections import Counter

def get_molecular_formula(name):
    """
    Returns the molecular formula for a given chemical name.
    This acts as a lookup table for the problem.
    """
    formulas = {
        # Reactants
        "1,1-dimethoxyethan-1-amine": "C4H11NO2",
        "but-3-en-2-ol": "C4H8O",
        "methanol": "CH4O",
        "(3R,4S)-3,4-dimethylhexa-1,5-diyne": "C8H10",
        "2-((vinyloxy)methyl)but-1-ene": "C7H12O",
        # Products from all options
        "6-methyl-3,4-dihydro-2H-pyran-2-amine": "C6H11NO",
        "(Z)-1-(but-2-en-2-yloxy)ethen-1-amine": "C6H11NO",
        "(1Z,2E)-1,2-diethylidenecyclobutane": "C8H12",
        "(3Z,4E)-3,4-diethylidenecyclobut-1-ene": "C8H10",
        "4-methylenehexanal": "C7H12O",
        "4-methylenehexan-1-ol": "C7H14O",
    }
    return formulas.get(name)

def parse_formula(formula_str):
    """
    Parses a molecular formula string (e.g., 'C4H11NO2') into a dictionary
    of atom counts (e.g., {'C': 4, 'H': 11, 'N': 1, 'O': 2}).
    """
    if not formula_str:
        return None
    # Use Counter for easier addition/subtraction
    return Counter(dict(re.findall(r'([A-Z][a-z]*)(\d*)', formula_str)))

def check_answer(option_string):
    """
    Checks the correctness of a given option string.
    The option string should be in the format "A = name1, B = name2, C = name3"
    """
    try:
        # Parse the product names from the option string
        products = dict(re.findall(r'([ABC])\s*=\s*([^,]+(?:,[^,]+)*)', option_string))
        product_A_name = products['A'].strip()
        product_B_name = products['B'].strip()
        product_C_name = products['C'].strip()
    except (KeyError, IndexError):
        return "Invalid option string format."

    # --- Check Reaction A ---
    reactant1_formula = parse_formula(get_molecular_formula("1,1-dimethoxyethan-1-amine"))
    reactant2_formula = parse_formula(get_molecular_formula("but-3-en-2-ol"))
    byproduct_formula = parse_formula(get_molecular_formula("methanol")) * 2
    expected_A_formula = reactant1_formula + reactant2_formula - byproduct_formula
    product_A_formula = parse_formula(get_molecular_formula(product_A_name))

    if expected_A_formula != product_A_formula:
        return f"Incorrect for Reaction A. Atom conservation failed. Expected formula {expected_A_formula}, but product '{product_A_name}' has formula {product_A_formula}."

    # --- Check Reaction B ---
    reactant_B_formula = parse_formula(get_molecular_formula("(3R,4S)-3,4-dimethylhexa-1,5-diyne"))
    product_B_formula = parse_formula(get_molecular_formula(product_B_name))

    if reactant_B_formula != product_B_formula:
        return f"Incorrect for Reaction B. Product must be an isomer of the reactant. Expected formula {reactant_B_formula}, but product '{product_B_name}' has formula {product_B_formula}."

    # --- Check Reaction C ---
    reactant_C_formula = parse_formula(get_molecular_formula("2-((vinyloxy)methyl)but-1-ene"))
    product_C_formula = parse_formula(get_molecular_formula(product_C_name))

    if reactant_C_formula != product_C_formula:
        return f"Incorrect for Reaction C. Product must be an isomer of the reactant. Expected formula {reactant_C_formula}, but product '{product_C_name}' has formula {product_C_formula}."

    # Check functional group for Claisen rearrangement
    if not (product_C_name.endswith("al") or product_C_name.endswith("one")):
        return f"Incorrect for Reaction C. A Claisen rearrangement should produce a carbonyl compound (aldehyde or ketone), but product '{product_C_name}' is not."

    return "Correct"

# The answer provided by the LLM corresponds to option D.
# Let's format it for the checker function.
llm_answer_string = "A = (Z)-1-(but-2-en-2-yloxy)ethen-1-amine, B = (3Z,4E)-3,4-diethylidenecyclobut-1-ene, C = 4-methylenehexanal"

# Run the check
result = check_answer(llm_answer_string)
print(result)
