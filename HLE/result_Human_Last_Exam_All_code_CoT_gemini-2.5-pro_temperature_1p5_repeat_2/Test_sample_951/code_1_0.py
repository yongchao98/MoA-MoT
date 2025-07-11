def solve_phonon_count():
    """
    Calculates the number of IR-active phonons for LiNiPO4 in different polarizations.
    """
    # 1. Define the irreducible representations for the D2h point group.
    irreps = ["Ag", "B1g", "B2g", "B3g", "Au", "B1u", "B2u", "B3u"]

    # 2. Define the character contributions for each Wyckoff site in LiNiPO4.
    # These decompositions are standard results from group theory for the Pnma space group.
    # gamma_4a is for 4 atoms on the 4a site (e.g., Li).
    gamma_4a = {"Ag": 1, "B1g": 2, "B2g": 1, "B3g": 2, "Au": 2, "B1u": 1, "B2u": 2, "B3u": 1}
    # gamma_4c is for 4 atoms on the 4c site (e.g., Ni, P, O1, O2).
    gamma_4c = {"Ag": 2, "B1g": 1, "B2g": 2, "B3g": 1, "Au": 1, "B1u": 2, "B2u": 1, "B3u": 2}
    # gamma_8d is for 8 atoms on the 8d site (e.g., O3).
    gamma_8d = {"Ag": 3, "B1g": 3, "B2g": 3, "B3g": 3, "Au": 3, "B1u": 3, "B2u": 3, "B3u": 3}

    # 3. Calculate the total representation for all 28 atoms (84 modes).
    # There is 1 species on 4a, 4 species on 4c, and 1 species on 8d.
    total_gamma = {irrep: 0 for irrep in irreps}
    for irrep in irreps:
        total_gamma[irrep] = gamma_4a[irrep] + (4 * gamma_4c[irrep]) + gamma_8d[irrep]

    # 4. Define the acoustic modes for the D2h point group.
    acoustic_modes = {"B1u": 1, "B2u": 1, "B3u": 1}

    # 5. Calculate the optical modes by subtracting the acoustic modes.
    optical_gamma = total_gamma.copy()
    for irrep, count in acoustic_modes.items():
        optical_gamma[irrep] -= count

    # 6. Determine the number of IR-active modes for each polarization.
    # E||x corresponds to B3u, E||y to B2u, E||z to B1u.
    num_x = optical_gamma["B3u"]
    num_y = optical_gamma["B2u"]
    num_z = optical_gamma["B1u"]

    # 7. Print the final result in the specified format.
    # The f-string directly includes the calculated numbers in the final output string.
    print(f"E||x: {num_x}, E||y: {num_y}, E||z: {num_z}")

solve_phonon_count()