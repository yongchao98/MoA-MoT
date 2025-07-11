import re
import collections

def calculate_6_311G_star_star_primitives():
    """
    Calculates the total number of primitive Gaussian functions
    for a given molecule with the 6-311G** basis set.
    """
    # --- You can change the molecular formula here ---
    # Examples: "H2O", "CH4", "C6H6", "NH3"
    molecular_formula = "CH3Cl"

    # Number of primitive Gaussians for different elements in the 6-311G** basis set.
    # This is derived from the standard definition of the basis set.
    # H: 5s (from 3,1,1 split) + 3p (from **) = 8
    # Period 2 (C,N,O..): 6s(core) + 20sp(valence) + 6d(**) = 32
    # Period 3 (Si,P,S..): Have a more complex core structure, summing to 48 primitives.
    primitives_per_atom = {
        # Period 1
        'H': 8,
        # Period 2
        'Li': 32, 'Be': 32, 'B': 32, 'C': 32, 'N': 32, 'O': 32, 'F': 32, 'Ne': 32,
        # Period 3
        'Na': 48, 'Mg': 48, 'Al': 48, 'Si': 48, 'P': 48, 'S': 48, 'Cl': 48, 'Ar': 48,
    }

    # Parse the molecular formula into a dictionary of atom counts.
    # e.g., "CH3Cl" -> {'C': 1, 'H': 3, 'Cl': 1}
    tokens = re.findall('([A-Z][a-z]?)(\d*)', molecular_formula)
    atom_counts = collections.defaultdict(int)
    for element, count in tokens:
        if not element in primitives_per_atom:
            print(f"Warning: Element '{element}' is not supported by this script. It will be ignored.")
            continue
        atom_counts[element] += int(count) if count else 1

    if not atom_counts:
        print("No supported atoms found in the formula.")
        return

    # Calculate the total number of primitives and create the equation string.
    total_primitives = 0
    equation_parts = []
    
    # Process atoms in a canonical order (C, H, then alphabetically)
    element_order = sorted(atom_counts.keys(), key=lambda x: ('CH'.find(x), x))

    for element in element_order:
        count = atom_counts[element]
        prims = primitives_per_atom[element]
        total_primitives += count * prims
        equation_parts.append(f"{count} * {prims} ({element})")

    # Print the final result
    print(f"Calculation for {molecular_formula} with 6-311G** basis set:")
    final_equation = " + ".join(equation_parts) + f" = {total_primitives}"
    print(final_equation)
    print(f"\nTotal primitive Gaussians = {total_primitives}")


# Run the calculation
calculate_6_311G_star_star_primitives()