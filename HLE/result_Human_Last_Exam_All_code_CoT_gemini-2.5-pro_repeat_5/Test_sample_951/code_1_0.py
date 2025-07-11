def solve_phonons():
    """
    Performs a factor group analysis for LiNiPO4 to determine the number of
    IR-active phonons for different polarizations.
    """

    # 1. Define crystal structure and group theory data
    
    # Point group irreps for D2h
    d2h_irreps = ["Ag", "B1g", "B2g", "B3g", "Au", "B1u", "B2u", "B3u"]

    # Atom counts on each Wyckoff site
    # 5 atom types (Li, Ni, P, O1, O2) are on 4c sites.
    # 1 atom type (O3) is on the 8d site.
    atoms_on_4c = 5
    atoms_on_8d = 1
    
    # Representation of the displacement vector (x,y,z) in the site symmetry group
    # For 4c site (Cs symmetry with mirror plane sigma_xz):
    # (x,z) are symmetric (A'), y is anti-symmetric (A'')
    disp_rep_cs = {"A'": 2, "A''": 1}
    # For 8d site (C1 symmetry): All directions are symmetric (A)
    disp_rep_c1 = {"A": 3}

    # Correlation tables from site symmetry group to D2h point group
    correlation_cs_xz_to_d2h = {
        "A'": ["Ag", "B2g", "B1u", "B3u"],
        "A''": ["B1g", "B3g", "Au", "B2u"]
    }
    correlation_c1_to_d2h = {
        "A": d2h_irreps
    }

    # 2. Calculate total vibrational modes (Gamma_total)
    
    total_modes = {irrep: 0 for irrep in d2h_irreps}

    # Add contributions from the 5 atom types on 4c sites
    for site_irrep, multiplicity in disp_rep_cs.items():
        for d2h_irrep in correlation_cs_xz_to_d2h[site_irrep]:
            total_modes[d2h_irrep] += atoms_on_4c * multiplicity
            
    # Add contributions from the 1 atom type on 8d sites
    for site_irrep, multiplicity in disp_rep_c1.items():
        for d2h_irrep in correlation_c1_to_d2h[site_irrep]:
            total_modes[d2h_irrep] += atoms_on_8d * multiplicity

    # 3. Subtract acoustic modes to get optical modes (Gamma_optic)
    
    optical_modes = total_modes.copy()
    acoustic_modes = ["B1u", "B2u", "B3u"]
    for mode in acoustic_modes:
        optical_modes[mode] -= 1

    # 4. Print the full decomposition of optical modes (the "equation")
    
    gamma_optic_parts = []
    for irrep in d2h_irreps:
        count = optical_modes[irrep]
        if count > 0:
            gamma_optic_parts.append(f"{count}{irrep}")
    
    print("Full decomposition of optical modes:")
    print("Γ_optic = " + " + ".join(gamma_optic_parts))
    print("-" * 30)

    # 5. Identify and count IR active modes for each polarization
    
    # In D2h: x -> B3u, y -> B2u, z -> B1u
    num_x = optical_modes["B3u"]
    num_y = optical_modes["B2u"]
    num_z = optical_modes["B1u"]

    # 6. Print the final result in the requested format
    
    print("Predicted number of IR-active phonons:")
    print(f"E||x: {num_x}, E||y: {num_y}, E||z: {num_z}")


solve_phonons()
<<<E||x: 12, E||y: 7, E||z: 12>>>