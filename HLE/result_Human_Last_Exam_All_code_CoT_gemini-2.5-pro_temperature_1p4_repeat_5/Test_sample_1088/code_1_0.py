import re
from collections import Counter

def calculate_primitives_6_311Gss(formula: str):
    """
    Calculates and explains the number of primitive Gaussians in a 6-311G**
    basis set for a given chemical formula.

    This function focuses on Hydrogen and 2nd-row atoms (Li-Ne), as the
    naming convention is most clearly defined for them.
    """
    # --- Define the rules for 6-311G** ---
    # 6:   The core orbitals are represented by 1 contracted Gaussian function (CGF)
    #      built from 6 primitive Gaussian functions (PGFs).
    # 311: The valence orbitals are split into three parts using 3, 1, and 1 PGFs.
    # **:  Polarization functions are added.
    #      - First *: Adds one set of d-functions to heavy atoms (6 primitives).
    #      - Second *: Adds one set of p-functions to Hydrogen atoms (3 primitives).

    # Rule for Hydrogen (1s valence)
    h_s_prims = 3 + 1 + 1
    h_p_prims = 3 * 1
    h_total = h_s_prims + h_p_prims

    # Rule for 2nd row heavy atoms (e.g., C) (1s core, 2s/2p valence)
    heavy_core_s = 6
    heavy_val_s = 3 + 1 + 1
    heavy_val_p = 3 * (3 + 1 + 1)
    heavy_pol_d = 6 * 1
    heavy_total = heavy_core_s + heavy_val_s + heavy_val_p + heavy_pol_d

    # --- Parse the chemical formula ---
    try:
        parsed_atoms = re.findall(r'([A-Z][a-z]?)(\d*)', formula)
        if not parsed_atoms or ''.join([s + c for s, c in parsed_atoms]) != formula:
            raise ValueError("Invalid chemical formula format.")

        atom_counts = Counter()
        for symbol, count_str in parsed_atoms:
            count = int(count_str) if count_str else 1
            atom_counts[symbol] += count
    except (ValueError, TypeError):
        print(f"Error: Could not parse the chemical formula '{formula}'.")
        return

    # --- Perform Calculation & Print Explanation ---
    total_primitives = 0
    final_equation_parts = []
    
    # Define which atoms are 2nd-row for this rule
    second_row_atoms = {'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne'}
    
    # Process elements in a consistent order (C, H, then others)
    sorted_elements = sorted(atom_counts.keys(), key=lambda x: (x not in ('C', 'H'), x))

    print(f"--- Analysis of Primitives in 6-311G** for {formula} ---")

    for element in sorted_elements:
        count = atom_counts[element]
        per_atom_total = 0
        
        if element == 'H':
            per_atom_total = h_total
            print(f"\nAtom Type: Hydrogen (H) - Found {count} atom(s)")
            print(f"  - Valence s-orbitals ('-311'): 3 + 1 + 1 = {h_s_prims} primitives")
            print(f"  - Polarization p-orbitals ('**'): 1 set * 3 functions = {h_p_prims} primitives")
            print(f"  - Total per H atom = {h_s_prims} + {h_p_prims} = {per_atom_total}")

        elif element in second_row_atoms:
            per_atom_total = heavy_total
            print(f"\nAtom Type: {element} (2nd Row) - Found {count} atom(s)")
            print(f"  - Core s-orbital ('6-'): {heavy_core_s} primitives")
            print(f"  - Valence s-orbital ('-311'): 3 + 1 + 1 = {heavy_val_s} primitives")
            print(f"  - Valence p-orbitals ('-311'): 3 orbitals * (3 + 1 + 1) = {heavy_val_p} primitives")
            print(f"  - Polarization d-orbitals ('**'): 1 set * 6 functions = {heavy_pol_d} primitives")
            print(f"  - Total per {element} atom = {heavy_core_s} + {heavy_val_s} + {heavy_val_p} + {heavy_pol_d} = {per_atom_total}")

        else:
            print(f"\nAtom Type: {element} - Found {count} atom(s)")
            print("  - NOTE: The simple 6-311G** counting rule is ambiguous for this element.")
            print("  - It will be excluded from the total.")
            
        if per_atom_total > 0:
            total_primitives += count * per_atom_total
            final_equation_parts.append(f"({count} * {per_atom_total})")

    # --- Print Final Equation ---
    if final_equation_parts:
        print("\n" + "="*55)
        print(f"Total primitives for {formula}:")
        equation_str = " + ".join(final_equation_parts)
        print(f"  Total = {equation_str} = {total_primitives}")
        print("="*55)

# Demonstrate the calculation for Methane (CH4), which contains both atom types.
calculate_primitives_6_311Gss('CH4')