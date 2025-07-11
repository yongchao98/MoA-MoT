import re
from collections import Counter

def calculate_primitive_gaussians(molecular_formula):
    """
    Calculates the total number of primitive Gaussian functions in a 6-311G**
    basis set for a given molecular formula.
    """
    # Number of primitive Gaussians per atom for the 6-311G** basis set.
    # H: (3+1+1)s + 1p(3) = 8
    # Period 2 (Li-Ne): 6(core) + (3+1+1)s + 3*(3+1+1)p + 1d(6) = 6 + 5 + 15 + 6 = 32
    # Period 3 (Na-Ar): 42 (from BSE) + 1d(6) = 48
    primitives_per_atom = {
        # Period 1
        'H': 8, 'He': 8,
        # Period 2
        'Li': 32, 'Be': 32, 'B': 32, 'C': 32, 'N': 32, 'O': 32, 'F': 32, 'Ne': 32,
        # Period 3
        'Na': 48, 'Mg': 48, 'Al': 48, 'Si': 48, 'P': 48, 'S': 48, 'Cl': 48, 'Ar': 48
    }

    # --- Parse the molecular formula ---
    # This regex finds atom symbols (e.g., H, C, Cl) followed by an optional number.
    parts = re.findall(r'([A-Z][a-z]?)(\d*)', molecular_formula)
    
    # Create a Counter dictionary to store the count of each atom.
    # e.g., 'C2H4O2' -> {'C': 2, 'H': 4, 'O': 2}
    atom_counts = Counter()
    for symbol, count in parts:
        count = int(count) if count else 1
        atom_counts[symbol] += count

    # --- Calculate and Print the Result ---
    total_primitives = 0
    calculation_steps = []
    
    print(f"Molecule: {molecular_formula}")
    print("Basis Set: 6-311G**")
    print("---")
    print("Calculation:")
    
    # Sort atoms for consistent output order (e.g., C, H, O)
    sorted_atoms = sorted(atom_counts.keys())

    for atom_symbol in sorted_atoms:
        if atom_symbol not in primitives_per_atom:
            print(f"Warning: Atom '{atom_symbol}' is not supported by this script. Skipping.")
            continue
            
        count = atom_counts[atom_symbol]
        primitives = primitives_per_atom[atom_symbol]
        sub_total = count * primitives
        total_primitives += sub_total
        
        step_str = f"{count} * {primitives}"
        calculation_steps.append(str(sub_total))
        
        print(f"{atom_symbol}: {count} atoms * {primitives} primitives/atom = {sub_total}")

    print("---")
    final_equation = " + ".join(calculation_steps)
    print(f"Total Primitive Gaussians = {final_equation} = {total_primitives}")


if __name__ == '__main__':
    # You can change this molecular formula to any molecule of your choice.
    # The parser can handle both condensed (CH3COOH) and simple (C2H4O2) formulas.
    molecular_formula = "C2H4O2"  # Acetic Acid
    calculate_primitive_gaussians(molecular_formula)