import re
from collections import Counter

def calculate_primitive_gaussians(formula: str):
    """
    Calculates and explains the total number of primitive Gaussian functions
    for a given molecule with the 6-311G** basis set.

    This script demonstrates the calculation for the molecule specified in the
    `formula_to_calculate` variable.
    """

    # --- Step 1: Define primitives per atom for 6-311G** ---
    # This is generally valid for 2nd-row elements (Li-Ne).
    # Heavy Atom (like C, N, O): 6(core) + 4*(3+1+1)(valence) + 6(d-pol) = 32
    # Hydrogen Atom: (3+1+1)(valence) + 3(p-pol) = 8
    primitives_per_atom = {
        'H': 8,
        'C': 32, 'N': 32, 'O': 32, 'F': 32, 'B': 32, 'Be': 32, 'Li': 32, 'Ne': 32
    }

    # --- Step 2: Parse the chemical formula to get atom counts ---
    try:
        atom_counts = Counter()
        # Find all atom-count pairs (e.g., 'C', 'H2', 'O')
        for atom, count_str in re.findall(r'([A-Z][a-z]*)(\d*)', formula):
            count = int(count_str) if count_str else 1
            atom_counts[atom] += count
    except Exception:
        print(f"Error: Could not parse the chemical formula '{formula}'.")
        print("Please use a standard format like 'H2O', 'CH4', or 'C6H6'.")
        return

    # --- Step 3: Calculate the total and prepare the explanation ---
    total_primitives = 0
    calculation_parts = []
    equation_parts = []
    value_parts = []

    # Sort atoms for consistent output (e.g., C, H, then others)
    sorted_atoms = sorted(atom_counts.keys(), key=lambda x: (x != 'C', x != 'H', x))

    for atom in sorted_atoms:
        count = atom_counts[atom]
        if atom not in primitives_per_atom:
            print(f"Warning: The calculation for atom '{atom}' is not defined.")
            continue

        primitives = primitives_per_atom[atom]
        total_primitives += count * primitives

        calculation_parts.append(f"({count} * {atom})")
        equation_parts.append(f"({count} * {primitives})")
        value_parts.append(str(count * primitives))

    # --- Step 4: Print the detailed explanation and result ---
    print(f"To find the number of primitive Gaussians in a 6-311G** basis set for {formula}:")

    print("\nFirst, we determine the number of primitives for each type of atom:")
    print("  - For a heavy atom (like C, N, O):")
    print("    - Core orbitals (1s) are represented by 6 primitives.")
    print("    - Valence orbitals (2s, 2p) are split into 3 parts (3+1+1), for 4 orbitals: 4 * (3 + 1 + 1) = 20 primitives.")
    print("    - The first '*' adds one set of d-polarization functions: 6 primitives.")
    print("    - Total per heavy atom = 6 + 20 + 6 = 32 primitives.")
    print("\n  - For a Hydrogen atom (H):")
    print("    - Valence orbital (1s) is split into 3 parts: 3 + 1 + 1 = 5 primitives.")
    print("    - The second '*' adds one set of p-polarization functions: 3 primitives.")
    print("    - Total per Hydrogen atom = 5 + 3 = 8 primitives.")

    print(f"\nNext, we apply this to the molecule {formula}, which has:")
    for atom in sorted_atoms:
        print(f"  - {atom_counts[atom]} atom(s) of {atom}")

    print("\nThe total number of primitive Gaussians is calculated as follows:")
    # Print the equation with atom symbols, e.g., Total = (1 * C) + (2 * H) + (1 * O)
    print(f"  Total = {' + '.join(calculation_parts)}")
    # Print the equation with primitive counts, e.g., Total = (1 * 32) + (2 * 8) + (1 * 32)
    print(f"        = {' + '.join(equation_parts)}")
    # Print the equation with intermediate results, e.g., Total = 32 + 16 + 32
    if len(value_parts) > 1:
        print(f"        = {' + '.join(value_parts)}")
    # Print the final result
    print(f"        = {total_primitives}")
    
    return total_primitives

if __name__ == '__main__':
    # We will use formaldehyde (CH2O) as our example molecule.
    formula_to_calculate = "CH2O"
    final_answer = calculate_primitive_gaussians(formula_to_calculate)
    print(f"\n<<<{final_answer}>>>")