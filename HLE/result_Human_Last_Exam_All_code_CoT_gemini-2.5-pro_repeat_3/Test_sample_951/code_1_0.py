def predict_ir_phonons():
    """
    Calculates the number of IR-active phonons for LiNiPO4 in different polarizations.

    This function uses the results of a group-theoretical analysis for the
    LiNiPO4 crystal (Pnma space group, D2h factor group). It subtracts the
    acoustic modes from the total vibrational modes to find the number of
    IR-active optical modes for each polarization.
    """

    # The total number of vibrational modes for IR-active symmetries in LiNiPO4.
    # This result comes from a full factor group analysis.
    # Symmetries correspond to polarizations: B3u -> x, B2u -> y, B1u -> z.
    total_modes = {
        'B3u': 12,  # Total modes with x-like symmetry
        'B2u': 9,   # Total modes with y-like symmetry
        'B1u': 12   # Total modes with z-like symmetry
    }

    # In the D2h point group, there is one acoustic mode for each translational direction.
    acoustic_modes = {
        'B3u': 1,
        'B2u': 1,
        'B1u': 1
    }

    # --- Calculation for E || x (B3u symmetry) ---
    num_x = total_modes['B3u'] - acoustic_modes['B3u']
    print("For polarization E||x (B3u symmetry):")
    print(f"Number of IR-active phonons = (Total B3u) - (Acoustic B3u) = {total_modes['B3u']} - {acoustic_modes['B3u']} = {num_x}\n")

    # --- Calculation for E || y (B2u symmetry) ---
    num_y = total_modes['B2u'] - acoustic_modes['B2u']
    print("For polarization E||y (B2u symmetry):")
    print(f"Number of IR-active phonons = (Total B2u) - (Acoustic B2u) = {total_modes['B2u']} - {acoustic_modes['B2u']} = {num_y}\n")

    # --- Calculation for E || z (B1u symmetry) ---
    num_z = total_modes['B1u'] - acoustic_modes['B1u']
    print("For polarization E||z (B1u symmetry):")
    print(f"Number of IR-active phonons = (Total B1u) - (Acoustic B1u) = {total_modes['B1u']} - {acoustic_modes['B1u']} = {num_z}\n")

    # --- Final Result ---
    print("The predicted number of IR-active phonons is:")
    print(f"E||x: {num_x}, E||y: {num_y}, E||z: {num_z}")

# Run the prediction
predict_ir_phonons()