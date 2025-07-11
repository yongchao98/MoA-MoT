import re

def calculate_molecular_weight(formula):
    """
    Calculates the molecular weight of a chemical formula and prints the detailed calculation.

    Args:
        formula (str): The chemical formula (e.g., 'C20H34O4Si').
    """
    atomic_weights = {
        'C': 12.011,
        'H': 1.008,
        'O': 15.999,
        'Si': 28.085,
        'N': 14.007,
        'P': 30.974,
        'S': 32.06,
        'Cl': 35.45,
        'F': 18.998,
        'Br': 79.904,
        'I': 126.90
    }

    # Use regex to find all element-count pairs in the formula
    # It handles elements with one or two letters and optional counts
    elements = re.findall(r'([A-Z][a-z]*)(\d*)', formula)

    total_mw = 0
    calculation_parts = []
    
    print(f"The reaction is an isomerization, so the product has the same molecular formula as the starting material.")
    print(f"Product Molecular Formula: {formula}\n")
    print("Molecular Weight Calculation:")

    for element, count in elements:
        # If count is not specified, it's 1
        count = int(count) if count else 1
        
        if element in atomic_weights:
            weight = atomic_weights[element]
            contribution = count * weight
            total_mw += contribution
            
            # Print contribution of each element
            print(f"Contribution from {element}: {count} * {weight} = {contribution:.3f}")
            calculation_parts.append(f"{contribution:.3f}")
        else:
            print(f"Warning: Atomic weight for element '{element}' not found.")
            return

    # Print the final equation
    equation_str = " + ".join(calculation_parts)
    print(f"\nTotal MW = {equation_str} = {total_mw:.3f}")

# The molecular formula of the starting material, and thus the product.
molecular_formula = 'C20H34O4Si'
calculate_molecular_weight(molecular_formula)