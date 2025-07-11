import re

def calculate_primitives_for_molecule():
    """
    Calculates the total number of primitive Gaussian functions in a 6-311G**
    basis set for a given molecule.
    """
    # The user did not specify a molecule, so we will use Caffeine (C8H10N4O2)
    # as a representative example. You can change the formula here.
    molecule_formula = "C8H10N4O2"

    print("The 6-311G** basis set defines the number of primitive Gaussians per atom as follows:")
    print("- For each Hydrogen (H) atom: 8 primitives")
    print("- For each second-row heavy atom (e.g., C, N, O): 32 primitives")
    print("-" * 50)

    # --- Formula Parsing ---
    # This regular expression finds element symbols (e.g., 'C', 'He') and their optional counts.
    tokens = re.findall(r'([A-Z][a-z]*)(\d*)', molecule_formula)
    atom_counts = {}
    for element, count_str in tokens:
        # If count is not specified, it defaults to 1.
        count = int(count_str) if count_str else 1
        # Add to existing count if element is repeated (e.g., in CH3COOH)
        atom_counts[element] = atom_counts.get(element, 0) + count

    # --- Primitive Calculation ---
    primitives_heavy = 32
    primitives_h = 8
    
    total_primitives = 0
    calculation_details = []
    final_equation_parts = []
    
    # Sort elements for a consistent output order (e.g., C, H, N, O)
    sorted_elements = sorted(atom_counts.keys())

    print(f"Calculating for the molecule: {molecule_formula}")
    print(f"Atom counts: {', '.join([f'{el}: {atom_counts[el]}' for el in sorted_elements])}\n")
    print("Breakdown of the calculation:")

    for element in sorted_elements:
        count = atom_counts[element]
        
        if element == 'H':
            primitives_per_atom = primitives_h
        else:
            # Assume any non-Hydrogen atom is a "heavy atom" for this basis set.
            primitives_per_atom = primitives_heavy
            
        subtotal = count * primitives_per_atom
        total_primitives += subtotal
        
        calculation_details.append(f"- {element}: {count} atoms * {primitives_per_atom} primitives/atom = {subtotal}")
        final_equation_parts.append(f"({count} * {primitives_per_atom})")

    for detail in calculation_details:
        print(detail)

    final_equation = " + ".join(final_equation_parts)
    print("\nFinal equation showing each number:")
    print(f"{final_equation} = {total_primitives}")


# Execute the function
calculate_primitives_for_molecule()