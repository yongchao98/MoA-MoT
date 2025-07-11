def calculate_ir_phonons():
    """
    Calculates the number of IR-active optical phonons for LiNiPO4
    based on its known group theoretical decomposition.
    """
    
    # The total decomposition of the 84 vibrational modes for LiNiPO4
    # into the irreducible representations of the D2h point group.
    # This is the result of a full Factor Group Analysis.
    # Format: { 'irrep_name': number_of_modes }
    total_modes = {
        'Ag': 12, 'B1g': 9, 'B2g': 12, 'B3g': 9,
        'Au': 9, 'B1u': 12, 'B2u': 9, 'B3u': 12
    }

    # Acoustic modes transform as x, y, z, which correspond to B3u, B2u, and B1u.
    acoustic_modes_irreps = ['B1u', 'B2u', 'B3u']

    # IR active modes have the same symmetry as the acoustic modes.
    # We calculate the number of *optical* IR-active modes by subtracting
    # the single acoustic mode from the total count for each relevant irrep.
    
    # For polarization E||x (corresponds to B3u symmetry)
    num_x = total_modes['B3u'] - 1
    
    # For polarization E||y (corresponds to B2u symmetry)
    num_y = total_modes['B2u'] - 1
    
    # For polarization E||z (corresponds to B1u symmetry)
    num_z = total_modes['B1u'] - 1
    
    # Print the result in the required format
    # Remember to show the final equation with each number explicitly.
    print(f"E||x: {num_x}, E||y: {num_y}, E||z: {num_z}")

if __name__ == "__main__":
    calculate_ir_phonons()
