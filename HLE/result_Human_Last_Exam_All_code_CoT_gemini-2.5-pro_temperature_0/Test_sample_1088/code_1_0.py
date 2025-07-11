def calculate_basis_set_primitives():
    """
    Calculates and explains the number of primitive Gaussians in the 6-311G**
    basis set for Hydrogen, a representative heavy atom (Carbon), and Methane (CH4).
    """
    print("The number of primitive Gaussians in a 6-311G** basis set depends on the atom.")
    print("Let's calculate it for Hydrogen and a representative heavy atom, Carbon.\n")

    # --- Calculation for Hydrogen (H) ---
    print("--- For a Hydrogen (H) atom ---")
    # The 1s orbital is valence, described by the '311' part.
    h_s_primitives = 3 + 1 + 1
    # The second '*' adds p-polarization functions (px, py, pz).
    h_p_pol_primitives = 3
    h_total = h_s_primitives + h_p_pol_primitives
    print("Valence 1s orbital primitives = 3 + 1 + 1 = 5")
    print("Polarization p-function primitives = 3")
    print(f"Total primitives for H = {h_s_primitives} + {h_p_pol_primitives} = {h_total}")
    print("-" * 35)

    # --- Calculation for a Heavy Atom (Carbon, C) ---
    print("\n--- For a second-row heavy atom (e.g., Carbon, C) ---")
    # The '6' describes the 1s core orbital.
    c_core_s_primitives = 6
    # The '311' describes the valence 2s and 2p orbitals.
    c_valence_s_primitives = 3 + 1 + 1
    # There are three 2p orbitals (px, py, pz).
    c_valence_p_primitives = 3 * (3 + 1 + 1)
    # The first '*' adds d-polarization functions (6 Cartesian d-functions).
    c_d_pol_primitives = 6
    c_total = c_core_s_primitives + c_valence_s_primitives + c_valence_p_primitives + c_d_pol_primitives
    print(f"Core 1s orbital primitives = {c_core_s_primitives}")
    print(f"Valence 2s orbital primitives = 3 + 1 + 1 = {c_valence_s_primitives}")
    print(f"Valence 2p orbitals primitives = 3 * (3 + 1 + 1) = {c_valence_p_primitives}")
    print(f"Polarization d-function primitives = {c_d_pol_primitives}")
    print(f"Total primitives for C = {c_core_s_primitives} + {c_valence_s_primitives} + {c_valence_p_primitives} + {c_d_pol_primitives} = {c_total}")
    print("-" * 35)

    # --- Example for Methane (CH4) ---
    print("\n--- Example: Total for Methane (CH4) ---")
    num_c = 1
    num_h = 4
    total_methane = (num_c * c_total) + (num_h * h_total)
    print(f"Methane has {num_c} Carbon atom and {num_h} Hydrogen atoms.")
    print(f"Total primitives = ({num_c} * primitives for C) + ({num_h} * primitives for H)")
    print(f"Total for CH4 = ({num_c} * {c_total}) + ({num_h} * {h_total}) = {total_methane}")


# Execute the calculation
calculate_basis_set_primitives()