def predict_ir_phonons():
    """
    Calculates the number of expected IR-active phonons for LiNiPO4.

    The calculation is based on the known factor group analysis for the
    olivine structure with space group Pnma (point group D2h).
    """

    # Step 1: Known decomposition of the total 84 vibrational modes for LiNiPO4
    # into the irreducible representations (irreps) of the D2h point group.
    total_modes_decomposition = {
        'Ag': 12, 'B1g': 9, 'B2g': 12, 'B3g': 9,
        'Au': 9, 'B1u': 12, 'B2u': 9, 'B3u': 12
    }

    # Step 2: Identify the irreps corresponding to IR activity for each polarization.
    # In the standard Pnma setting, the electric field vector E transforms as:
    # E || x  -> B3u symmetry
    # E || y  -> B2u symmetry
    # E || z  -> B1u symmetry
    ir_active_irreps = {
        'x': 'B3u',
        'y': 'B2u',
        'z': 'B1u'
    }

    # Step 3: Calculate the number of IR-active *optical* phonons.
    # For each IR-active symmetry, one mode is acoustic (translation). The rest are optical.
    # We subtract 1 from the total number of modes for each of these symmetries.

    ir_x_irrep = ir_active_irreps['x']
    num_b3u_total = total_modes_decomposition[ir_x_irrep]
    num_ir_x = num_b3u_total - 1

    ir_y_irrep = ir_active_irreps['y']
    num_b2u_total = total_modes_decomposition[ir_y_irrep]
    num_ir_y = num_b2u_total - 1

    ir_z_irrep = ir_active_irreps['z']
    num_b1u_total = total_modes_decomposition[ir_z_irrep]
    num_ir_z = num_b1u_total - 1

    # Step 4: Print the results, showing the calculation for each polarization.
    print("Calculation of IR-active optical phonons for LiNiPO4:")
    print(f"For E||x (Symmetry {ir_x_irrep}): Total Modes({num_b3u_total}) - Acoustic Modes(1) = {num_ir_x}")
    print(f"For E||y (Symmetry {ir_y_irrep}): Total Modes({num_b2u_total}) - Acoustic Modes(1) = {num_ir_y}")
    print(f"For E||z (Symmetry {ir_z_irrep}): Total Modes({num_b1u_total}) - Acoustic Modes(1) = {num_ir_z}")
    
    print("\n---")
    print("Predicted number of phonons in polarized IR spectra:")
    print(f"E||x: {num_ir_x}, E||y: {num_ir_y}, E||z: {num_ir_z}")

# Run the prediction
predict_ir_phonons()