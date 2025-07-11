def calculate_primitives_for_6_311G_star_star():
    """
    Calculates and explains the number of primitive Gaussians in a 6-311G**
    basis set for Hydrogen, Carbon, and a Methane molecule.
    """

    print("--- Calculation for a Hydrogen (H) atom ---")
    # For H, the 1s orbital is valence.
    # Valence 's' primitives from '311': 3 + 1 + 1 = 5
    h_s_primitives = 3 + 1 + 1
    # Polarization from the second '*' adds one set of p-functions (px, py, pz).
    # Each is a single uncontracted primitive. So, 3 * 1 = 3
    h_p_primitives = 3 * 1
    total_h = h_s_primitives + h_p_primitives
    print(f"A Hydrogen atom has {h_s_primitives} (s-type) + {h_p_primitives} (p-type) primitive Gaussians.")
    print(f"Equation: {h_s_primitives} + {h_p_primitives} = {total_h}\n")


    print("--- Calculation for a Carbon (C) atom ---")
    # Carbon is a "heavy atom".
    # Core '1s' primitives from '6-': 6
    c_core_s_primitives = 6
    # Valence '2s' primitives from '-311': 3 + 1 + 1 = 5
    c_valence_s_primitives = 3 + 1 + 1
    # Valence '2p' primitives from '-311': 3 sets (px, py, pz) of (3 + 1 + 1) = 3 * 5 = 15
    c_valence_p_primitives = 3 * (3 + 1 + 1)
    # Polarization from the first '*' adds one set of d-functions.
    # A standard set has 6 functions (dxy, dxz, dyz, d_x2-y2, d_z2, etc.),
    # each being a single primitive. So, 6 * 1 = 6
    c_d_primitives = 6 * 1
    total_c = c_core_s_primitives + c_valence_s_primitives + c_valence_p_primitives + c_d_primitives
    print(f"A Carbon atom has {c_core_s_primitives} (core s) + {c_valence_s_primitives} (valence s) + {c_valence_p_primitives} (valence p) + {c_d_primitives} (d-type) primitives.")
    print(f"Equation: {c_core_s_primitives} + {c_valence_s_primitives} + {c_valence_p_primitives} + {c_d_primitives} = {total_c}\n")

    print("--- Example: Total for a Methane (CH4) molecule ---")
    # Methane has 1 Carbon atom and 4 Hydrogen atoms.
    total_methane = total_c + 4 * total_h
    print(f"Methane (CH4) has 1 Carbon atom and 4 Hydrogen atoms.")
    print(f"Total primitives = (Primitives for C) + 4 * (Primitives for H)")
    print(f"Final Equation: {total_c} + 4 * {total_h} = {total_methane}")

calculate_primitives_for_6_311G_star_star()