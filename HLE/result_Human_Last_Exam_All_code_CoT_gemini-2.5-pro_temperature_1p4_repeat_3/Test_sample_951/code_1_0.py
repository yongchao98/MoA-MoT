def calculate_ir_phonons():
    """
    Calculates the number of IR-active phonons for LiNiPO4 based on
    its known factor group analysis.
    """

    # Total number of modes for each irreducible representation (irrep) for LiNiPO4 (Pnma)
    total_modes = {
        'B1u': 12,  # Corresponds to polarization E||z in the standard Pnma setting
        'B2u': 9,   # Corresponds to polarization E||y
        'B3u': 12,  # Corresponds to polarization E||x
    }

    # There is one acoustic mode for each of the IR-active irreps.
    acoustic_modes = {
        'B1u': 1,
        'B2u': 1,
        'B3u': 1,
    }

    # --- Calculation for E||x (B3u symmetry) ---
    total_b3u = total_modes['B3u']
    acoustic_b3u = acoustic_modes['B3u']
    ir_b3u = total_b3u - acoustic_b3u
    print(f"Calculation for E||x (B3u modes):")
    print(f"Total B3u modes = {total_b3u}")
    print(f"Acoustic B3u modes = {acoustic_b3u}")
    print(f"Number of IR active phonons = {total_b3u} - {acoustic_b3u} = {ir_b3u}\n")
    
    num_x = ir_b3u

    # --- Calculation for E||y (B2u symmetry) ---
    total_b2u = total_modes['B2u']
    acoustic_b2u = acoustic_modes['B2u']
    ir_b2u = total_b2u - acoustic_b2u
    print(f"Calculation for E||y (B2u modes):")
    print(f"Total B2u modes = {total_b2u}")
    print(f"Acoustic B2u modes = {acoustic_b2u}")
    print(f"Number of IR active phonons = {total_b2u} - {acoustic_b2u} = {ir_b2u}\n")
    
    num_y = ir_b2u

    # --- Calculation for E||z (B1u symmetry) ---
    total_b1u = total_modes['B1u']
    acoustic_b1u = acoustic_modes['B1u']
    ir_b1u = total_b1u - acoustic_b1u
    print(f"Calculation for E||z (B1u modes):")
    print(f"Total B1u modes = {total_b1u}")
    print(f"Acoustic B1u modes = {acoustic_b1u}")
    print(f"Number of IR active phonons = {total_b1u} - {acoustic_b1u} = {ir_b1u}\n")

    num_z = ir_b1u

    # Print the final result in the requested format
    print("Final result:")
    print(f"E||x: {num_x}, E||y: {num_y}, E||z: {num_z}")


if __name__ == '__main__':
    calculate_ir_phonons()
