def calculate_6_311G_star_star_primitives():
    """
    Calculates and explains the number of primitive Gaussians in a 6-311G** basis set
    for Hydrogen, a representative heavy atom (Carbon), and Methane (CH4).
    """

    # --- Calculation for a Hydrogen Atom ---
    print("1. For a Hydrogen atom (H):")
    # Valence 1s orbital (311 split)
    h_s_primitives = 3 + 1 + 1
    # Polarization p-functions (second '*') -> px, py, pz
    h_p_pol_primitives = 3
    # Total for H
    h_total_primitives = h_s_primitives + h_p_pol_primitives
    print(f"   - The 1s orbital is described by a '311' split: 3 + 1 + 1 = {h_s_primitives} primitives.")
    print(f"   - The second '*' adds one set of p-polarization functions: {h_p_pol_primitives} primitives.")
    print(f"   - Total for Hydrogen = {h_s_primitives} (s) + {h_p_pol_primitives} (p) = {h_total_primitives} primitive Gaussians.\n")

    # --- Calculation for a 2nd Row Heavy Atom (e.g., Carbon) ---
    print("2. For a second-row heavy atom (e.g., Carbon, C):")
    # Core 1s orbital ('6' split)
    heavy_core_s_primitives = 6
    # Valence 2s orbital ('311' split)
    heavy_valence_s_primitives = 3 + 1 + 1
    # Valence 2p orbitals (px, py, pz), each with a '311' split
    heavy_valence_p_primitives = 3 * (3 + 1 + 1)
    # Polarization d-functions (first '*') -> 6 Cartesian d-functions
    heavy_d_pol_primitives = 6
    # Total for C
    heavy_total_primitives = heavy_core_s_primitives + heavy_valence_s_primitives + heavy_valence_p_primitives + heavy_d_pol_primitives
    print(f"   - The 1s core orbital is described by '6': {heavy_core_s_primitives} primitives.")
    print(f"   - The 2s valence orbital is described by '311': 3 + 1 + 1 = {heavy_valence_s_primitives} primitives.")
    print(f"   - The three 2p valence orbitals are each described by '311': 3 * (3 + 1 + 1) = {heavy_valence_p_primitives} primitives.")
    print(f"   - The first '*' adds one set of d-polarization functions: {heavy_d_pol_primitives} primitives.")
    print(f"   - Total for Carbon = {heavy_core_s_primitives} (1s) + {heavy_valence_s_primitives} (2s) + {heavy_valence_p_primitives} (2p) + {heavy_d_pol_primitives} (d) = {heavy_total_primitives} primitive Gaussians.\n")

    # --- Example Calculation for Methane (CH4) ---
    print("3. Example: For a Methane molecule (CH4):")
    # 1 Carbon atom
    num_c = 1
    # 4 Hydrogen atoms
    num_h = 4
    total_ch4_primitives = (num_c * heavy_total_primitives) + (num_h * h_total_primitives)
    print(f"   - Methane has {num_c} Carbon atom and {num_h} Hydrogen atoms.")
    print(f"   - Total Primitives = ({num_c} * {heavy_total_primitives}) + ({num_h} * {h_total_primitives}) = {total_ch4_primitives} primitive Gaussians.")
    
    return total_ch4_primitives

if __name__ == '__main__':
    # The final answer required is a single number. We provide it for a common molecule, CH4.
    final_answer = calculate_6_311G_star_star_primitives()
    # print(f"\nFinal Answer for Methane (CH4): {final_answer}")
    # The final answer is wrapped according to the instruction.
    print(f"\n<<<{final_answer}>>>")
