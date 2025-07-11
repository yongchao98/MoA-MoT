import re
from collections import Counter

def parse_formula(formula):
    """Parses a chemical formula string into a dictionary of element counts."""
    return Counter({m[0]: int(m[1] or 1) for m in re.finditer(r'([A-Z][a-z]*)(\d*)', formula)})

def format_formula(atom_counts):
    """Formats a dictionary of element counts into a chemical formula string."""
    # Order elements C, H, then alphabetically for others
    elements = sorted(atom_counts.keys(), key=lambda x: ('Z' if x not in 'CH' else x))
    if 'H' in elements:
        elements.insert(0, elements.pop(elements.index('H')))
    if 'C' in elements:
        elements.insert(0, elements.pop(elements.index('C')))
    
    formula_str = ""
    for el in elements:
        count = atom_counts[el]
        formula_str += el
        if count > 1:
            formula_str += str(count)
    return formula_str

def verify_reaction_stoichiometry(reaction_name, start_formula, alkyne_formula, loss_formula, product_formula_given):
    """
    Verifies the molecular formula of the product by simulating the reaction steps.
    The mechanism involves formation of a münchnone (loss of H2O), cycloaddition,
    and then extrusion of another molecule (e.g., CO2).
    """
    print(f"--- Verifying Stoichiometry for {reaction_name} ---")
    
    start_counts = parse_formula(start_formula)
    h2o_counts = parse_formula("H2O")
    alkyne_counts = parse_formula(alkyne_formula)
    loss_counts = parse_formula(loss_formula)
    
    # Step 1: Form münchnone intermediate by losing H2O from the starting carboxylic acid
    munchone_counts = start_counts - h2o_counts
    munchone_formula = format_formula(munchone_counts)
    print(f"Starting Material ({start_formula}) - H2O => Münchnone Intermediate ({munchone_formula})")

    # Step 2: [3+2] Cycloaddition with the alkyne
    cycloadduct_counts = munchone_counts + alkyne_counts
    cycloadduct_formula = format_formula(cycloadduct_counts)
    print(f"Münchnone ({munchone_formula}) + Alkyne ({alkyne_formula}) => Primary Cycloadduct ({cycloadduct_formula})")
    
    # Step 3: Extrusion of CO2 to form the final product
    final_product_counts = cycloadduct_counts - loss_counts
    final_product_formula_calc = format_formula(final_product_counts)
    print(f"Primary Cycloadduct ({cycloadduct_formula}) - {loss_formula} => Calculated Final Product ({final_product_formula_calc})")
    
    # Step 4: Compare the calculated formula with the given HRMS data
    print(f"\nCalculated Formula: {final_product_formula_calc}")
    print(f"Given Formula:      {product_formula_given}")
    if final_product_formula_calc == product_formula_given:
        print("Result: The calculated formula MATCHES the provided HRMS data.\n")
    else:
        print("Result: The calculated formula DOES NOT MATCH the provided HRMS data.\n")

# --- Define molecular formulas from the problem ---
# Reaction 1/2 starting material is N-acetyl-N-methyl-alanine
start_rxn1_2 = "C6H11NO3" 
product_A_given = "C9H13NO2"

# Reaction 3 starting material is N-acetylproline
start_rxn3 = "C7H11NO3"
product_B_given = "C10H13NO2"

# Common reagent (methyl propiolate) and lost molecule (carbon dioxide)
methyl_propiolate = "C4H4O2"
co2_loss = "CO2"

# --- Run the verifications ---
verify_reaction_stoichiometry("Reaction 1 & 2 / Product A", start_rxn1_2, methyl_propiolate, co2_loss, product_A_given)
verify_reaction_stoichiometry("Reaction 3 / Product B", start_rxn3, methyl_propiolate, co2_loss, product_B_given)