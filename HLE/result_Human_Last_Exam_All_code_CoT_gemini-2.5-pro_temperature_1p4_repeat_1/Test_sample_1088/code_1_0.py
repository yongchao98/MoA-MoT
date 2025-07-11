import re
from collections import Counter

def calculate_primitives(formula: str):
    """
    Calculates the total number of primitive Gaussian functions for a given
    molecule with the 6-311G** basis set.
    """
    # The number of primitive Gaussians for each element in the 6-311G** basis set.
    # The name is shorthand:
    # - For H and He, it's 6-311G**.
    # - For all other atoms (Li-Ar), it's 6-311G*.
    # Data is sourced from the Basis Set Exchange (bse.pnl.gov).
    primitive_counts = {
        'H': 8, 'He': 8,
        'Li': 17, 'Be': 17,
        'B': 32, 'C': 32, 'N': 32, 'O': 32, 'F': 32, 'Ne': 32,
        'Na': 42, 'Mg': 42,
        'Al': 51, 'Si': 51, 'P': 51, 'S': 51, 'Cl': 51, 'Ar': 51,
    }

    # Use a regular expression to find all atoms and their optional counts
    atoms_found = re.findall(r'([A-Z][a-z]*)(\d*)', formula)
    
    if not atoms_found:
        print(f"Error: Could not parse the formula '{formula}'. Please use a standard format (e.g., CH4, H2O, C6H12O6).")
        return

    # Create a Counter dictionary to handle multiple occurrences of the same element
    atom_counts = Counter()
    for symbol, count_str in atoms_found:
        if symbol not in primitive_counts:
            print(f"Error: Basis set data for element '{symbol}' is not available in this script.")
            return
        count = int(count_str) if count_str else 1
        atom_counts[symbol] += count

    total_primitives = 0
    equation_parts = []
    calculation_parts = []

    # Sort atoms alphabetically for a consistent output order (e.g., C, H, O)
    for symbol in sorted(atom_counts.keys()):
        count = atom_counts[symbol]
        prims_per_atom = primitive_counts[symbol]
        
        term_total = count * prims_per_atom
        total_primitives += term_total
        
        # Build the parts for the final equation string
        equation_parts.append(f"({count} * {prims_per_atom})")
        calculation_parts.append(str(term_total))

    print(f"Calculation for {formula} with the 6-311G** basis set:")
    print("-" * 40)
    
    # Print the breakdown for each element
    for i, symbol in enumerate(sorted(atom_counts.keys())):
        print(f"For {symbol}: {atom_counts[symbol]} atom(s) * {primitive_counts[symbol]} primitives/atom = {calculation_parts[i]}")

    # Print the final combined equation showing all numbers
    final_equation_str = " + ".join(equation_parts)
    final_calc_str = " + ".join(calculation_parts)
    
    print("-" * 40)
    print(f"Total = {final_equation_str}")
    print(f"      = {final_calc_str}")
    print(f"      = {total_primitives}")
    print("-" * 40)


if __name__ == '__main__':
    # Since no molecule was specified, we use methane (CH4) as an example.
    # You can change this to any other formula.
    molecule_formula = "CH4"
    calculate_primitives(molecule_formula)
    # Another example: Water (H2O)
    # calculate_primitives("H2O")
    # Another example: Chlorine Trifluoride (ClF3)
    # calculate_primitives("ClF3")

<<<64>>>