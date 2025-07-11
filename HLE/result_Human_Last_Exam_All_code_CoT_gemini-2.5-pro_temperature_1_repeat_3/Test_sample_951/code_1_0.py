def solve_phonon_modes():
    """
    Calculates the number of IR-active phonons for LiNiPO4 in different polarizations.
    """
    # The symmetry operations of the D2h point group in a standard order
    ops = ['E', 'C2z', 'C2y', 'C2x', 'i', 'sigma_z', 'sigma_y', 'sigma_x']

    # D2h character table. The rows correspond to irreducible representations (irreps).
    # The columns correspond to the symmetry operations in the 'ops' list.
    D2h_char_table = {
        'Ag':  [1,  1,  1,  1,  1,  1,  1,  1],
        'B1g': [1,  1, -1, -1,  1,  1, -1, -1],
        'B2g': [1, -1,  1, -1,  1, -1,  1, -1],
        'B3g': [1, -1, -1,  1,  1, -1, -1,  1],
        'Au':  [1,  1,  1,  1, -1, -1, -1, -1],
        'B1u': [1,  1, -1, -1, -1, -1,  1,  1],  # Transforms as T_z
        'B2u': [1, -1,  1, -1, -1,  1, -1,  1],  # Transforms as T_y
        'B3u': [1, -1, -1,  1, -1,  1,  1, -1],  # Transforms as T_x
    }

    # Number of atoms left unmoved by each symmetry operation for LiNiPO4 (Pnma)
    # Total atoms N = 28.
    unmoved_atoms = {
        'E': 28, 'C2z': 0, 'C2y': 0, 'C2x': 0,  # Screw axes move all atoms
        'i': 4,                                 # 4 Li atoms on inversion sites (4a)
        'sigma_z': 0, 'sigma_x': 0,             # Glide planes move all atoms
        'sigma_y': 16                           # 16 atoms on mirror planes (4c)
    }

    # Character of the vector representation (how a vector transforms)
    vec_char = {
        'E': 3, 'C2z': -1, 'C2y': -1, 'C2x': -1,
        'i': -3, 'sigma_z': 1, 'sigma_y': 1, 'sigma_x': 1
    }

    # Calculate the character of the total mechanical representation
    mech_char_list = [unmoved_atoms[op] * vec_char[op] for op in ops]

    # Decompose the mechanical representation into irreps
    group_order = 8
    total_modes = {}
    for irrep, char_list in D2h_char_table.items():
        dot_product = sum(i * j for i, j in zip(char_list, mech_char_list))
        total_modes[irrep] = dot_product // group_order
    
    # Subtract acoustic modes (B1u + B2u + B3u) to get optical modes
    optical_modes = total_modes.copy()
    optical_modes['B1u'] -= 1
    optical_modes['B2u'] -= 1
    optical_modes['B3u'] -= 1

    # Format the optical mode decomposition string
    optical_mode_str_parts = []
    for irrep in ['Ag', 'B1g', 'B2g', 'B3g', 'Au', 'B1u', 'B2u', 'B3u']:
        count = optical_modes[irrep]
        if count > 0:
            optical_mode_str_parts.append(f"{count}{irrep}")
    
    print("Decomposition of optical phonon modes:")
    print(f"Gamma_opt = {' + '.join(optical_mode_str_parts)}")
    print("-" * 30)

    # Assign the number of modes to each polarization
    # Assuming standard orientation: x || a, y || b, z || c
    # E||x -> B3u, E||y -> B2u, E||z -> B1u
    num_x = optical_modes['B3u']
    num_y = optical_modes['B2u']
    num_z = optical_modes['B1u']

    print("Predicted number of IR-active phonons:")
    print(f"E||x: {num_x}, E||y: {num_y}, E||z: {num_z}")

solve_phonon_modes()