import re

def calculate_6_311Gss_primitives(molecule_formula: str):
    """
    Calculates the total number of primitive Gaussian functions for a given molecule
    using the 6-311G** basis set.

    This function supports H and common 2nd-row elements (C, N, O, F).
    """
    # Primitives per atom for the 6-311G** basis set
    # H: 1s(3,1,1) + p(pol) = (3+1+1) + 3 = 8
    # 2nd Row (e.g., C,N,O): 1s(core, 6) + 2s(val, 3,1,1) + 2p(val, 3,1,1) + d(pol)
    #                      = 6 + (3+1+1) + 3*(3+1+1) + 6 = 6 + 5 + 15 + 6 = 32
    primitives_per_atom = {
        'H': 8,
        # Add other 2nd-row elements, they have the same count
        'C': 32,
        'N': 32,
        'O': 32,
        'F': 32,
    }

    # Find all element symbols followed by an optional number
    atom_counts = re.findall(r'([A-Z][a-z]*)(\d*)', molecule_formula)

    total_primitives = 0
    calculation_parts = []
    equation_parts = []
    
    print(f"Analyzing formula: {molecule_formula}")

    for element, count_str in atom_counts:
        if element not in primitives_per_atom:
            print(f"Warning: Atom '{element}' is not supported in this script. Skipping.")
            continue

        # If count is not specified, it's 1
        count = int(count_str) if count_str else 1
        
        primitives = primitives_per_atom[element]
        subtotal = count * primitives
        total_primitives += subtotal

        calculation_parts.append(f"({count} {element} atoms * {primitives} primitives/{element})")
        equation_parts.append(f"({count} * {primitives})")

    if not calculation_parts:
        print("Could not find any supported atoms in the formula.")
        return

    # Print the detailed breakdown of the calculation
    print("Calculation: " + " + ".join(calculation_parts))
    
    # Print the equation with numbers
    equation_str = " + ".join(equation_parts)
    subtotal_values = [str(eval(part)) for part in equation_parts] # e.g., ['16', '32']
    
    # Final equation format: (2 * 8) + (1 * 32) = 16 + 32 = 48
    print(f"Equation: {equation_str} = {' + '.join(subtotal_values)} = {total_primitives}")
    
    print(f"\nTotal number of primitive Gaussians in 6-311G** for {molecule_formula} is: {total_primitives}")
    return total_primitives

if __name__ == '__main__':
    # --- User can change the molecule formula here ---
    # Example molecules: "CH4", "C6H6", "H2O2", "NH3"
    formula = "H2O"
    
    calculate_6_311Gss_primitives(formula)