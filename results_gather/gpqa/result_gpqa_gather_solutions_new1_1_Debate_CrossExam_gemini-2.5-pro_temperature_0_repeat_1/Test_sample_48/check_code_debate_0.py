def get_molecular_formula(name):
    """
    Returns the molecular formula (C, H, N, O atoms) for a given compound name.
    This is a simplified lookup for the specific compounds in this problem.
    """
    formulas = {
        # Reactants
        "(3R,4S)-3,4-dimethylhexa-1,5-diyne": (8, 10, 0, 0),
        "2-((vinyloxy)methyl)but-1-ene": (7, 12, 0, 1),

        # Products for B
        "(3Z,4E)-3,4-diethylidenecyclobut-1-ene": (8, 10, 0, 0),
        "(1Z,2E)-1,2-diethylidenecyclobutane": (8, 12, 0, 0),

        # Products for C
        "4-methylenehexanal": (7, 12, 0, 1),
        "4-methylenehexan-1-ol": (7, 14, 0, 1),
    }
    return formulas.get(name)

def get_functional_group_type(name):
    """Identifies the functional group type from the IUPAC name suffix."""
    if name.endswith("al") or name.endswith("one"):
        return "carbonyl"
    if name.endswith("ol"):
        return "alcohol"
    return "other"

# Define the four multiple-choice options
options = {
    "A": {
        "B": "(3Z,4E)-3,4-diethylidenecyclobut-1-ene",
        "C": "4-methylenehexanal"
    },
    "B": {
        "B": "(1Z,2E)-1,2-diethylidenecyclobutane",
        "C": "4-methylenehexanal"
    },
    "C": {
        "B": "(3Z,4E)-3,4-diethylidenecyclobut-1-ene",
        "C": "4-methylenehexan-1-ol"
    },
    "D": {
        "B": "(1Z,2E)-1,2-diethylidenecyclobutane",
        "C": "4-methylenehexan-1-ol"
    }
}

# The final answer provided by the LLM to be checked
llm_answer = "A"

# --- Verification Logic ---

# Check Constraint 1: Reaction B must be an isomerization
reactant_b_formula = get_molecular_formula("(3R,4S)-3,4-dimethylhexa-1,5-diyne")
survivors_from_b_check = set()
for option_key, products in options.items():
    product_b_formula = get_molecular_formula(products["B"])
    if product_b_formula == reactant_b_formula:
        survivors_from_b_check.add(option_key)

# Check Constraint 2: Reaction C must be an isomerization to a carbonyl
reactant_c_formula = get_molecular_formula("2-((vinyloxy)methyl)but-1-ene")
survivors_from_c_check = set()
for option_key, products in options.items():
    product_c_formula = get_molecular_formula(products["C"])
    product_c_group = get_functional_group_type(products["C"])
    
    if product_c_formula == reactant_c_formula and product_c_group == "carbonyl":
        survivors_from_c_check.add(option_key)

# Determine the final correct option by finding the intersection
final_survivors = survivors_from_b_check.intersection(survivors_from_c_check)

# Check if the LLM's answer is correct
if len(final_survivors) == 1 and llm_answer in final_survivors:
    print("Correct")
else:
    error_messages = []
    # Explain why the LLM's answer is wrong, if it is
    if llm_answer not in survivors_from_b_check:
        reactant_formula_str = f"C{reactant_b_formula[0]}H{reactant_b_formula[1]}"
        product_b = options[llm_answer]["B"]
        product_formula = get_molecular_formula(product_b)
        product_formula_str = f"C{product_formula[0]}H{product_formula[1]}"
        error_messages.append(f"Reaction B constraint is not satisfied: The product '{product_b}' (formula {product_formula_str}) is not an isomer of the reactant (formula {reactant_formula_str}).")

    if llm_answer not in survivors_from_c_check:
        reactant_formula_str = f"C{reactant_c_formula[0]}H{reactant_c_formula[1]}O"
        product_c = options[llm_answer]["C"]
        product_formula = get_molecular_formula(product_c)
        product_formula_str = f"C{product_formula[0]}H{product_formula[1]}O"
        product_group = get_functional_group_type(product_c)
        
        if product_formula != reactant_c_formula:
             error_messages.append(f"Reaction C constraint is not satisfied: The product '{product_c}' (formula {product_formula_str}) is not an isomer of the reactant (formula {reactant_formula_str}).")
        if product_group != "carbonyl":
             error_messages.append(f"Reaction C constraint is not satisfied: The product '{product_c}' is an {product_group}, not a carbonyl compound as expected from a Claisen rearrangement.")
    
    if not error_messages:
         error_messages.append(f"The provided answer '{llm_answer}' is incorrect because the constraints point to {final_survivors} as the only valid option(s).")

    print("Incorrect. " + " ".join(error_messages))
