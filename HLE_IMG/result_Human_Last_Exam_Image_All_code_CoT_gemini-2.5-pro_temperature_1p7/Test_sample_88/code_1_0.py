import re

def calculate_molecular_weight(formula):
    """
    Calculates the molecular weight of a chemical compound given its formula.
    
    Args:
        formula (str): The chemical formula string (e.g., 'C9H14N2O2').
        
    Returns:
        float: The calculated molecular weight.
    """
    atomic_masses = {
        'C': 12.011,
        'H': 1.008,
        'N': 14.007,
        'O': 15.999,
    }
    
    # Regex to find elements (e.g., C, H, N, O) and their counts
    pattern = re.compile('([A-Z][a-z]*)(\d*)')
    elements = pattern.findall(formula)
    
    total_mw = 0.0
    for (element_symbol, count_str) in elements:
        # If count is an empty string, it defaults to 1
        count = int(count_str) if count_str else 1
        
        if element_symbol in atomic_masses:
            total_mw += atomic_masses[element_symbol] * count
        else:
            raise ValueError(f"Element {element_symbol} not found in atomic mass data.")
            
    return total_mw

# Define the molecular formulas based on the proposed pathway for Product B
sm_formula = "C9H14N2O2"       # Starting Material
mp_formula = "C4H4O2"        # Methyl Propiolate
meoh_formula = "CH4O"         # Methanol (eliminated)
product_b_formula = "C12H14N2O3"  # Product B

# Calculate the molecular weights of all species
sm_mw = calculate_molecular_weight(sm_formula)
mp_mw = calculate_molecular_weight(mp_formula)
meoh_mw = calculate_molecular_weight(meoh_formula)
product_b_mw = calculate_molecular_weight(product_b_formula)

# Verify the pathway by checking mass conservation
pathway_result = sm_mw + mp_mw - meoh_mw

# Print the final equation with all the numbers
print("Verifying the pathway to Product B by calculating molecular weights:")
print(f"Proposed Reaction: Starting Material ({sm_formula}) + Methyl Propiolate ({mp_formula}) - Methanol ({meoh_formula}) -> Product B ({product_b_formula})")
print("\nEquation with calculated molecular weights:")
print(f"{sm_mw:.3f} + {mp_mw:.3f} - {meoh_mw:.3f} = {pathway_result:.3f}")
print(f"\nThe calculated molecular weight of Product B from its formula is: {product_b_mw:.3f}")

# The final answer is the molecular weight of product B.
print(f"\nFinal Answer (Molecular Weight of Product B): {product_b_mw:.3f}")
