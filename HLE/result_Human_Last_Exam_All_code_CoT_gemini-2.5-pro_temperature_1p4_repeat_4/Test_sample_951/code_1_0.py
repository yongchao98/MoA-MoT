def calculate_ir_phonons():
    """
    Calculates and prints the number of IR-active phonons for LiNiPO4.

    The calculation is based on the factor group analysis for the Pnma space group.
    """

    # 1. Total mechanical representation for LiNiPO4 (Pnma, Z=4).
    # This decomposition of the 84 total modes is taken from established literature.
    total_modes = {
        'Ag': 12, 'B1g': 9, 'B2g': 12, 'B3g': 9,
        'Au': 9,  'B1u': 12, 'B2u': 9,  'B3u': 12
    }

    # 2. Acoustic modes in the D2h point group.
    # B3u -> x, B2u -> y, B1u -> z
    acoustic_modes = {'B1u': 1, 'B2u': 1, 'B3u': 1}

    # 3. Map polarization to irreducible representations.
    # By convention for orthorhombic crystals: x || a, y || b, z || c.
    ir_map = {
        'E||x': 'B3u',
        'E||y': 'B2u',
        'E||z': 'B1u'
    }

    print("Calculation for the number of IR-active optical phonons:")
    
    # E||x (B3u)
    repre_x = ir_map['E||x']
    num_x_total = total_modes[repre_x]
    num_x_acoustic = acoustic_modes.get(repre_x, 0)
    num_x_optical = num_x_total - num_x_acoustic
    print(f"For E||x ({repre_x} modes): {num_x_total} (total) - {num_x_acoustic} (acoustic) = {num_x_optical}")

    # E||y (B2u)
    repre_y = ir_map['E||y']
    num_y_total = total_modes[repre_y]
    num_y_acoustic = acoustic_modes.get(repre_y, 0)
    num_y_optical = num_y_total - num_y_acoustic
    print(f"For E||y ({repre_y} modes): {num_y_total} (total) - {num_y_acoustic} (acoustic) = {num_y_optical}")

    # E||z (B1u)
    repre_z = ir_map['E||z']
    num_z_total = total_modes[repre_z]
    num_z_acoustic = acoustic_modes.get(repre_z, 0)
    num_z_optical = num_z_total - num_z_acoustic
    print(f"For E||z ({repre_z} modes): {num_z_total} (total) - {num_z_acoustic} (acoustic) = {num_z_optical}")

    print("\nPredicted number of phonons in polarized IR spectra:")
    # Print the final result in the requested format.
    print(f"E||x: {num_x_optical}, E||y: {num_y_optical}, E||z: {num_z_optical}")

if __name__ == "__main__":
    calculate_ir_phonons()