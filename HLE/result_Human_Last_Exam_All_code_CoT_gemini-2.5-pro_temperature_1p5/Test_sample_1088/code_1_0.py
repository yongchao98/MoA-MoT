import re

def calculate_total_primitives(formula):
    """
    Calculates the total number of primitive Gaussian functions for a given
    molecule using the 6-311G** basis set.
    """
    # A dictionary mapping atomic symbols to their PGF count in a 6-311G** basis set.
    # Counts are based on standard definitions used in computational chemistry packages.
    pgf_counts = {
        # Row 1
        'H': 8, 'He': 5,
        # Row 2
        'Li': 17, 'Be': 17,
        'B': 32, 'C': 32, 'N': 32, 'O': 32, 'F': 32, 'Ne': 32,
        # Row 3 (use more complex, accurate counts)
        'Na': 45, 'Mg': 45, 'Al': 45, 'Si': 45, 'P': 45, 'S': 45, 'Cl': 45, 'Ar': 45,
    }

    # Regex to parse chemical formulas e.g., H2O, C6H12O6
    atom_pattern = re.compile(r'([A-Z][a-z]?)(\d*)')
    
    # Parse the formula string and consolidate atom counts
    try:
        atom_list = atom_pattern.findall(formula)
        if not atom_list or ''.join([symbol + (count or '1') for symbol, count in atom_list]) != formula:
            raise ValueError("Invalid chemical formula format.")
        
        parsed_atoms = {}
        for symbol, count_str in atom_list:
            if symbol not in pgf_counts:
                print(f"Error: Atom '{symbol}' is not supported in this script.")
                return
            count = int(count_str) if count_str else 1
            parsed_atoms[symbol] = parsed_atoms.get(symbol, 0) + count
    except (TypeError, ValueError) as e:
        print(f"Could not parse formula '{formula}': {e}")
        return

    total_primitives = 0
    equation_parts = []
    
    # Sort atoms alphabetically for consistent output order
    for symbol in sorted(parsed_atoms.keys()):
        count = parsed_atoms[symbol]
        primitives_per_atom = pgf_counts[symbol]
        total_primitives += count * primitives_per_atom
        equation_parts.append(f"{count} * {primitives_per_atom}")
        
    print(f"Calculation for 6-311G** primitives in {formula}:")
    print(f"The total is the sum of (atom count * primitives per atom).")
    # This line prints the final equation with all numbers.
    print(f"Equation: {' + '.join(equation_parts)} = {total_primitives}")
    print(f"\nTotal number of primitive Gaussians for {formula} is {total_primitives}.")
    
    return total_primitives

# --- Main execution ---
# We will calculate the number of primitives for a water molecule (H2O) as an example.
# The result for just H2O will be returned as the final answer.
final_answer = calculate_total_primitives("H2O")

# The final answer will be printed in the special format below.
# print(f"<<<{final_answer}>>>")