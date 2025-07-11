def calculate_ir_phonons():
    """
    Calculates the number of IR-active phonons for different polarizations
    in an olivine LiNiPO4 crystal (Pnma space group).
    """

    # Step 1 & 2: Define the total vibrational modes from factor group analysis
    # for the olivine structure (Pnma, D2h point group).
    # The crystal has 28 atoms in the primitive cell, giving 3*28 = 84 total modes.
    # The decomposition into irreducible representations (irreps) is known:
    total_modes = {
        'Ag': 12, 'B1g': 9, 'B2g': 12, 'B3g': 9,
        'Au': 9, 'B1u': 12, 'B2u': 9, 'B3u': 12
    }
    print("Step 1: Total vibrational modes at the Gamma point for Pnma olivine structure:")
    print(f"Γ_total = {total_modes}\n")

    # Step 3: Identify and subtract the 3 acoustic modes.
    # In D2h, acoustic modes transform as B1u(z) + B2u(y) + B3u(x).
    acoustic_modes = {'B1u': 1, 'B2u': 1, 'B3u': 1}
    print("Step 2: Identify the acoustic modes to be subtracted:")
    print(f"Γ_acoustic = {acoustic_modes}\n")

    # Calculate optical modes by subtracting acoustic modes from total modes.
    optical_modes = total_modes.copy()
    for irrep, count in acoustic_modes.items():
        optical_modes[irrep] -= count

    print("Step 3: Calculate the optical modes (Γ_optical = Γ_total - Γ_acoustic):")
    print(f"Γ_optical = {optical_modes}\n")

    # Step 4: Identify which modes are IR-active and their polarization.
    # This is determined by the D2h character table.
    # B1u is polarized along z, B2u along y, and B3u along x.
    # 'g' (gerade) and 'Au' modes are IR-inactive.
    ir_polarization_map = {
        'B3u': 'x',  # Corresponds to E || x
        'B2u': 'y',  # Corresponds to E || y
        'B1u': 'z'   # Corresponds to E || z
    }
    print("Step 4: Identify IR-active modes and their polarization:")
    print(f"Modes polarized along x-axis (E||x): B3u")
    print(f"Modes polarized along y-axis (E||y): B2u")
    print(f"Modes polarized along z-axis (E||z): B1u\n")
    
    # Step 5: Count the number of phonons for each polarization.
    num_x = optical_modes.get('B3u', 0)
    num_y = optical_modes.get('B2u', 0)
    num_z = optical_modes.get('B1u', 0)

    # Final Result
    print("Final result:")
    # The final print statement constructs the output string as requested.
    print(f"E||x: {num_x}, E||y: {num_y}, E||z: {num_z}")


calculate_ir_phonons()