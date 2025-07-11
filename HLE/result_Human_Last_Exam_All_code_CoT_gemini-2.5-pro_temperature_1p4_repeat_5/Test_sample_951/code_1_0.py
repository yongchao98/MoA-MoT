def predict_ir_phonons():
    """
    Calculates the number of expected IR-active phonons for olivine LiNiPO4.

    This calculation is based on the results of a factor group analysis for the
    Pnma (D2h) space group with the specific atomic occupancies of LiNiPO4.
    """

    print("--- Analysis of IR-Active Phonons in LiNiPO4 ---")
    print("Crystal Structure: Olivine")
    print("Space Group: Pnma (No. 62)")
    print("Factor Group: D2h\n")

    # From factor group analysis for the LiNiPO4 structure (28 atoms/cell),
    # the total number of modes belonging to the IR-active irreducible
    # representations (irreps) are known.
    # The irrep B3u corresponds to polarization E||x (or a-axis).
    # The irrep B2u corresponds to polarization E||y (or b-axis).
    # The irrep B1u corresponds to polarization E||z (or c-axis).
    total_modes = {
        'B3u': 12,  # For x-polarization
        'B2u': 9,   # For y-polarization
        'B1u': 11   # For z-polarization -- wait, my manual calc was 12. Let me recheck.
    }
    
    # Re-checking the sum from my thought process:
    # B1u: 1 (from 4a) + 4*2 (from 4c) + 3 (from 8d) = 1 + 8 + 3 = 12
    # B2u: 2 (from 4a) + 4*1 (from 4c) + 3 (from 8d) = 2 + 4 + 3 = 9
    # B3u: 1 (from 4a) + 4*2 (from 4c) + 3 (from 8d) = 1 + 8 + 3 = 12
    # So the total modes should be 12, 9, 12. Let's correct the dict.
    
    total_modes = {
        'B3u': 12,  # For x-polarization
        'B2u': 9,   # For y-polarization
        'B1u': 12   # For z-polarization
    }


    # Each IR-active irrep includes one acoustic mode corresponding to
    # the translation of the whole crystal. These are not optical phonons.
    acoustic_modes = {
        'B3u': 1,
        'B2u': 1,
        'B1u': 1
    }

    # The number of observable IR phonons is the number of optical modes,
    # found by subtracting the acoustic modes from the total modes.
    num_phonons_x = total_modes['B3u'] - acoustic_modes['B3u']
    num_phonons_y = total_modes['B2u'] - acoustic_modes['B2u']
    num_phonons_z = total_modes['B1u'] - acoustic_modes['B1u']
    
    print("Calculation for E||x (B3u symmetry):")
    print(f"Total B3u modes = {total_modes['B3u']}")
    print(f"Acoustic B3u modes = {acoustic_modes['B3u']}")
    print(f"Number of optical phonons = {total_modes['B3u']} - {acoustic_modes['B3u']} = {num_phonons_x}\n")

    print("Calculation for E||y (B2u symmetry):")
    print(f"Total B2u modes = {total_modes['B2u']}")
    print(f"Acoustic B2u modes = {acoustic_modes['B2u']}")
    print(f"Number of optical phonons = {total_modes['B2u']} - {acoustic_modes['B2u']} = {num_phonons_y}\n")
    
    print("Calculation for E||z (B1u symmetry):")
    print(f"Total B1u modes = {total_modes['B1u']}")
    print(f"Acoustic B1u modes = {acoustic_modes['B1u']}")
    print(f"Number of optical phonons = {total_modes['B1u']} - {acoustic_modes['B1u']} = {num_phonons_z}\n")

    # Format the final answer as requested
    final_answer = f"E||x: {num_phonons_x}, E||y: {num_phonons_y}, E||z: {num_phonons_z}"
    print("--- Final Prediction ---")
    print(final_answer)

# Execute the function to get the prediction
predict_ir_phonons()