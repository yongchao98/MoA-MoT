def calculate_6_311G_star_star_primitives():
    """
    Calculates and prints the number of primitive Gaussians in a 6-311G**
    basis set for a methane (CH4) molecule.
    """
    # --- Calculation for a Heavy Atom (Carbon) ---
    c_core = 6
    c_valence_s = 3 + 1 + 1
    c_valence_p = 3 * (3 + 1 + 1)
    c_polarization_d = 6
    total_c = c_core + c_valence_s + c_valence_p + c_polarization_d
    
    print("--- Calculation for Carbon (C) in 6-311G** ---")
    print(f"Core 1s primitives ('6-'): {c_core}")
    print(f"Valence 2s primitives ('-311'): 3 + 1 + 1 = {c_valence_s}")
    print(f"Valence 2p primitives ('-311'): 3 * (3 + 1 + 1) = {c_valence_p}")
    print(f"Polarization d-functions ('*'): {c_polarization_d}")
    print(f"Total for Carbon = {c_core} + {c_valence_s} + {c_valence_p} + {c_polarization_d} = {total_c}")
    print("-" * 20)

    # --- Calculation for a Hydrogen Atom ---
    h_valence_s = 3 + 1 + 1
    h_polarization_p = 3
    total_h = h_valence_s + h_polarization_p
    
    print("--- Calculation for Hydrogen (H) in 6-311G** ---")
    print(f"Valence 1s primitives ('-311'): 3 + 1 + 1 = {h_valence_s}")
    print(f"Polarization p-functions ('**'): {h_polarization_p}")
    print(f"Total for Hydrogen = {h_valence_s} + {h_polarization_p} = {total_h}")
    print("-" * 20)

    # --- Total Calculation for Methane (CH4) ---
    num_c = 1
    num_h = 4
    total_ch4 = (num_c * total_c) + (num_h * total_h)

    print("--- Total for Methane (CH4) ---")
    print(f"Equation: ({num_c} * Primitives_C) + ({num_h} * Primitives_H)")
    print(f"Total primitives = ({num_c} * {total_c}) + ({num_h} * {total_h}) = {total_ch4}")

# Execute the calculation
calculate_6_311G_star_star_primitives()

# The final answer for a CH4 molecule is 64.