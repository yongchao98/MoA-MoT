def calculate_primitives_6_311g_star_star():
    """
    Calculates the number of primitive Gaussian functions for a molecule
    with the 6-311G** basis set.

    This function uses methane (CH4) as a representative molecule.
    """

    # --- Define the components of the 6-311G** basis set ---

    # For any atom's valence s and p orbitals
    valence_split = [3, 1, 1]
    valence_sum = sum(valence_split)

    # --- Carbon (heavy atom) ---
    num_c_atoms = 1
    # Core (1s orbital)
    c_core_primitives = 6
    # Valence s-orbital (2s)
    c_valence_s_primitives = valence_sum
    # Valence p-orbitals (2px, 2py, 2pz)
    c_valence_p_primitives = 3 * valence_sum
    # Polarization d-functions (from the first '*')
    # A standard set of d-functions has 6 components (d_xx, d_yy, d_zz, d_xy, d_xz, d_yz)
    c_polarization_primitives = 6

    # Total for one Carbon atom
    c_total = c_core_primitives + c_valence_s_primitives + c_valence_p_primitives + c_polarization_primitives

    # --- Hydrogen ---
    num_h_atoms = 4
    # Hydrogen has no core electrons; its 1s orbital is treated as valence.
    # Valence s-orbital (1s)
    h_valence_s_primitives = valence_sum
    # Polarization p-functions (from the second '*')
    # A set of p-functions has 3 components (px, py, pz)
    h_polarization_primitives = 3

    # Total for one Hydrogen atom
    h_total = h_valence_s_primitives + h_polarization_primitives

    # --- Methane (CH4) Molecule ---
    total_primitives = (num_c_atoms * c_total) + (num_h_atoms * h_total)

    # --- Print the explanation and results ---
    print("Calculating the number of primitive Gaussians in a 6-311G** basis set for a methane (CH4) molecule.")
    print("\n--- For one Carbon atom (a heavy atom) ---")
    print(f"- Core (1s orbital): represented by {c_core_primitives} primitives.")
    print(f"- Valence s (2s orbital): represented by {valence_split[0]} + {valence_split[1]} + {valence_split[2]} = {c_valence_s_primitives} primitives.")
    print(f"- Valence p (2p orbitals): 3 orbitals * ({valence_split[0]} + {valence_split[1]} + {valence_split[2]}) = {c_valence_p_primitives} primitives.")
    print(f"- Polarization (d-functions): represented by {c_polarization_primitives} primitives.")
    print(f"Equation: {c_core_primitives} + {c_valence_s_primitives} + {c_valence_p_primitives} + {c_polarization_primitives} = {c_total}")
    print(f"Total for one Carbon atom = {c_total} primitives.")


    print("\n--- For one Hydrogen atom ---")
    print(f"- Valence (1s orbital): represented by {valence_split[0]} + {valence_split[1]} + {valence_split[2]} = {h_valence_s_primitives} primitives.")
    print(f"- Polarization (p-functions): represented by {h_polarization_primitives} primitives.")
    print(f"Equation: {h_valence_s_primitives} + {h_polarization_primitives} = {h_total}")
    print(f"Total for one Hydrogen atom = {h_total} primitives.")

    print("\n--- For one Methane (CH4) molecule ---")
    print(f"Total primitives = ({num_c_atoms} * primitives for Carbon) + ({num_h_atoms} * primitives for Hydrogen)")
    print(f"Equation: {num_c_atoms} * ({c_total}) + {num_h_atoms} * ({h_total}) = {total_primitives}")
    print(f"\nThe total number of primitive Gaussians for CH4 with the 6-311G** basis set is {total_primitives}.")


if __name__ == "__main__":
    calculate_primitives_6_311g_star_star()