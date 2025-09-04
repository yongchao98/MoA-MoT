import re
from collections import Counter

def get_molecular_formula(name):
    """
    Returns the molecular formula for a given chemical name.
    In a real-world scenario, this would involve a complex cheminformatics library.
    For this problem, we can use a pre-populated dictionary based on the structures.
    """
    formulas = {
        # --- Reactants ---
        "1,1-dimethoxyethan-1-amine": "C4H11NO2",
        "but-3-en-2-ol": "C4H8O",
        "(3R,4S)-3,4-dimethylhexa-1,5-diyne": "C8H10",
        "2-((vinyloxy)methyl)but-1-ene": "C7H12O",
        # --- Products for A ---
        "6-methyl-3,4-dihydro-2H-pyran-2-amine": "C6H11NO",
        "(Z)-1-(but-2-en-2-yloxy)ethen-1-amine": "C6H11NO",
        # --- Products for B ---
        "(1Z,2E)-1,2-diethylidenecyclobutane": "C8H12",
        "(3Z,4E)-3,4-diethylidenecyclobut-1-ene": "C8H10",
        # --- Products for C ---
        "4-methylenehexanal": "C7H12O",
        "4-methylenehexan-1-ol": "C7H14O",
    }
    return formulas.get(name)

def parse_formula(formula_str):
    """Parses a chemical formula string into a Counter of atom counts."""
    if not formula_str:
        return Counter()
    return Counter({
        atom: int(count) if count else 1
        for atom, count in re.findall(r'([A-Z][a-z]*)(\d*)', formula_str)
    })

def check_reaction_A(product_A_name):
    """
    Checks if the stoichiometry for reaction A is plausible.
    The reaction is a condensation, likely eliminating two methanol (CH4O) molecules.
    """
    reactant1_formula = parse_formula(get_molecular_formula("1,1-dimethoxyethan-1-amine"))
    reactant2_formula = parse_formula(get_molecular_formula("but-3-en-2-ol"))
    product_A_formula = parse_formula(get_molecular_formula(product_A_name))
    
    total_reactants = reactant1_formula + reactant2_formula
    
    # The expected byproduct is two molecules of methanol (2 * CH4O = C2H8O2)
    byproduct_formula = parse_formula("C2H8O2")
    
    expected_product_formula = total_reactants - byproduct_formula
    
    if product_A_formula == expected_product_formula:
        return True, ""
    else:
        return False, f"Stoichiometry for product A ('{product_A_name}') is incorrect. Expected {expected_product_formula} but got {product_A_formula} after accounting for byproducts."

def check_reaction_B(product_B_name):
    """
    Checks if reaction B is a valid isomerization.
    The product must have the same molecular formula as the reactant.
    """
    reactant_formula = parse_formula(get_molecular_formula("(3R,4S)-3,4-dimethylhexa-1,5-diyne"))
    product_B_formula = parse_formula(get_molecular_formula(product_B_name))
    
    if reactant_formula == product_B_formula:
        return True, ""
    else:
        return False, f"Stoichiometry for product B ('{product_B_name}') is incorrect. It is not an isomer of the reactant. Reactant is {reactant_formula}, product is {product_B_formula}."

def check_reaction_C(product_C_name):
    """
    Checks if reaction C (Claisen rearrangement) is valid.
    1. It must be an isomerization (same molecular formula).
    2. The product must be a carbonyl compound (aldehyde or ketone).
    """
    # Check 1: Isomerization
    reactant_formula = parse_formula(get_molecular_formula("2-((vinyloxy)methyl)but-1-ene"))
    product_C_formula = parse_formula(get_molecular_formula(product_C_name))
    
    if reactant_formula != product_C_formula:
        return False, f"Stoichiometry for product C ('{product_C_name}') is incorrect. It is not an isomer of the reactant. Reactant is {reactant_formula}, product is {product_C_formula}."
        
    # Check 2: Functional Group
    # A Claisen rearrangement of an allyl vinyl ether yields a gamma,delta-unsaturated carbonyl.
    if not (product_C_name.endswith("al") or product_C_name.endswith("one")):
        return False, f"Functional group for product C ('{product_C_name}') is incorrect. A Claisen rearrangement should produce a carbonyl (aldehyde/ketone), not an alcohol or other group."
        
    return True, ""

def check_answer():
    """
    Checks the provided answer D against the chemical constraints.
    """
    # The final answer from the LLM is D.
    llm_answer = {
        "A": "(Z)-1-(but-2-en-2-yloxy)ethen-1-amine",
        "B": "(3Z,4E)-3,4-diethylidenecyclobut-1-ene",
        "C": "4-methylenehexanal"
    }

    # Check product A
    is_A_correct, reason_A = check_reaction_A(llm_answer["A"])
    if not is_A_correct:
        return f"Incorrect. Reason: {reason_A}"

    # Check product B
    is_B_correct, reason_B = check_reaction_B(llm_answer["B"])
    if not is_B_correct:
        return f"Incorrect. Reason: {reason_B}"

    # Check product C
    is_C_correct, reason_C = check_reaction_C(llm_answer["C"])
    if not is_C_correct:
        return f"Incorrect. Reason: {reason_C}"

    # If all checks pass
    return "Correct"

# Run the check and print the result
result = check_answer()
print(result)