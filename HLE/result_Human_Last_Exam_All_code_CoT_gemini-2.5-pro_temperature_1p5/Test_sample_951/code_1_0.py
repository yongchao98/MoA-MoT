def calculate_ir_phonons_linipo4():
    """
    Calculates and explains the number of expected IR-active phonons
    for olivine LiNiPO4 based on group theory.
    """
    print("Group Theory Analysis for Olivine LiNiPO4:")
    print("------------------------------------------")
    print("Material: LiNiPO4")
    print("Space Group: Pnma (No. 62), Point Group: D2h")
    print("Atoms per unit cell (N): 28")
    print("Total vibrational modes (3N): 84")
    print("\nThe decomposition of total modes into irreducible representations (irreps) is:")
    print("Γ_total = 12Ag + 9B1g + 12B2g + 9B3g + 9Au + 12B1u + 9B2u + 12B3u\n")


    # From factor group analysis, the total number of modes for each IR-active irrep
    total_modes_by_irrep = {
        'B1u': 12,  # Corresponds to z-polarization
        'B2u': 9,   # Corresponds to y-polarization
        'B3u': 12   # Corresponds to x-polarization
    }

    # There is 1 acoustic mode for each IR-active irrep
    acoustic_modes_count = 1

    print("Calculating IR-active optical phonons for each polarization:")
    print("---------------------------------------------------------")
    print("Note: IR-active optical modes = Total modes of a given symmetry - Acoustic modes of that symmetry.\n")

    # Polarization E || x (a-axis)
    irrep_x = 'B3u'
    total_x = total_modes_by_irrep[irrep_x]
    optical_x = total_x - acoustic_modes_count
    print(f"For E||x (a-axis), activity corresponds to {irrep_x} modes.")
    print(f"Calculation for E||x: {total_x} - {acoustic_modes_count} = {optical_x}\n")

    # Polarization E || y (b-axis)
    irrep_y = 'B2u'
    total_y = total_modes_by_irrep[irrep_y]
    optical_y = total_y - acoustic_modes_count
    print(f"For E||y (b-axis), activity corresponds to {irrep_y} modes.")
    print(f"Calculation for E||y: {total_y} - {acoustic_modes_count} = {optical_y}\n")

    # Polarization E || z (c-axis)
    irrep_z = 'B1u'
    total_z = total_modes_by_irrep[irrep_z]
    optical_z = total_z - acoustic_modes_count
    print(f"For E||z (c-axis), activity corresponds to {irrep_z} modes.")
    print(f"Calculation for E||z: {total_z} - {acoustic_modes_count} = {optical_z}\n")

    print("Final Answer:")
    print(f"E||x: {optical_x}, E||y: {optical_y}, E||z: {optical_z}")

# Execute the calculation
calculate_ir_phonons_linipo4()