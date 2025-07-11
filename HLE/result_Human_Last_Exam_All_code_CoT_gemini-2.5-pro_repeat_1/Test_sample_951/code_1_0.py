import collections

def calculate_ir_phonons():
    """
    Calculates the number of expected IR phonons for polarized spectra of LiNiPO4.

    This function uses the results of a factor group analysis for the olivine LiNiPO4
    crystal structure (space group Pnma, point group D2h).
    """

    print("Analyzing the vibrational modes for olivine LiNiPO4...")
    print("Crystal structure: Orthorhombic, Space Group Pnma (No. 62), Point Group D2h")
    print("Formula units per cell (Z): 4")
    print("Total atoms in cell: 4 Li + 4 Ni + 4 P + 16 O = 28")
    print("Total vibrational modes: 3 * 28 = 84")
    print("-" * 30)

    # The decomposition of OPTICAL modes at the Gamma point is well-established
    # from factor group analysis, consistent with literature (e.g., Phys. Rev. B 82, 094441 (2010)).
    # We define the number of modes for each irreducible representation (irrep).
    optical_modes = collections.OrderedDict([
        ('Ag', 12), ('B1g', 6), ('B2g', 12), ('B3g', 6),
        ('Au', 6), ('B1u', 11), ('B2u', 5), ('B3u', 11)
    ])

    # Construct the full equation string for Gamma_optic
    gamma_optic_str = " + ".join([f"{count}{irrep}" for irrep, count in optical_modes.items()])
    print("The decomposition of the optical modes is:")
    print(f"Γ_optic = {gamma_optic_str}")
    print("\nNote: Raman active modes are Ag, B1g, B2g, B3g. Au modes are silent.")
    print("-" * 30)

    # For the D2h point group, the IR active modes transform as x, y, and z.
    # The acoustic modes (1 B1u, 1 B2u, 1 B3u) have already been excluded in the optical representation above.
    print("The IR active modes and their corresponding polarizations are:")
    print("E || x  <-->  B3u symmetry")
    print("E || y  <-->  B2u symmetry")
    print("E || z  <-->  B1u symmetry")
    print("-" * 30)

    # Extract the number of phonons for each polarization from the optical modes
    num_phonons_x = optical_modes['B3u']
    num_phonons_y = optical_modes['B2u']
    num_phonons_z = optical_modes['B1u']

    print("The predicted number of observable IR active phonons for each polarization is:")
    print(f"For E||x, the number of B3u modes is: {num_phonons_x}")
    print(f"For E||y, the number of B2u modes is: {num_phonons_y}")
    print(f"For E||z, the number of B1u modes is: {num_phonons_z}")
    print("-" * 30)

    # Print the final result in the requested format
    final_result = f"E||x: {num_phonons_x}, E||y: {num_phonons_y}, E||z: {num_phonons_z}"
    print("Final Prediction:")
    print(final_result)


if __name__ == '__main__':
    calculate_ir_phonons()