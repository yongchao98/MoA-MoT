import re

def calculate_primitive_gaussians(formula: str):
    """
    Calculates the total number of primitive Gaussian functions for a given
    molecular formula using the 6-311G** basis set.

    Args:
        formula: A string representing the molecular formula (e.g., "C9H8O4").
    """
    # Number of primitive Gaussians for each element in the 6-311G** basis set.
    # Values are derived from standard quantum chemistry library (pyscf) definitions.
    pgf_counts = {
        # Period 1
        'H': 6, 'He': 5,
        # Period 2
        'Li': 17, 'Be': 17, 'B': 17, 'C': 17, 'N': 17, 'O': 17, 'F': 17, 'Ne': 17,
        # Period 3
        'Na': 21, 'Mg': 21, 'Al': 21, 'Si': 21, 'P': 21, 'S': 21, 'Cl': 21, 'Ar': 21,
    }

    # Use regex to find all atom-count pairs in the formula
    # e.g., "C9H8O4" -> [('C', '9'), ('H', '8'), ('O', '4')]
    atom_list = re.findall('([A-Z][a-z]*)(\d*)', formula)

    if not atom_list:
        print(f"Error: Could not parse the formula '{formula}'.")
        print("Please use a valid format, e.g., H2O, C6H6.")
        return

    total_primitives = 0
    calculation_terms = []
    equation_terms = []

    print(f"Calculation for {formula} with 6-311G** basis set:")
    for symbol, count_str in atom_list:
        # Default to 1 if no number follows the atom symbol (e.g., 'O' in 'H2O')
        count = int(count_str) if count_str else 1

        if symbol not in pgf_counts:
            print(f"Error: Atom '{symbol}' is not supported in this script.")
            return

        primitives_per_atom = pgf_counts[symbol]
        term_total = count * primitives_per_atom
        total_primitives += term_total
        
        calculation_terms.append(f"For {symbol}: {count} * {primitives_per_atom} = {term_total}")
        equation_terms.append(str(term_total))

    # Print the breakdown for each element
    for term in calculation_terms:
        print(term)

    # Print the final summation equation
    equation_str = " + ".join(equation_terms)
    print(f"Total = {equation_str} = {total_primitives}")
    
    # Return the final answer in the specified format
    print(f"\n<<<{total_primitives}>>>")


# --- Main execution ---
# You can change this formula to calculate for a different molecule.
molecular_formula = "C9H8O4"  # Aspirin
calculate_primitive_gaussians(molecular_formula)
