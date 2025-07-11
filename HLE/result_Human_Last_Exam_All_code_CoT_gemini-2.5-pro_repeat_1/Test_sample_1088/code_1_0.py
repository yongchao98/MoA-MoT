def calculate_and_print_primitives():
    """
    Calculates and explains the number of primitive Gaussians in a 6-311G**
    basis set for representative atoms and an example molecule (Methane).
    """
    print("The number of primitive Gaussians in a 6-311G** basis set depends on the atom.")
    print("Here is the breakdown for common atoms:\n")

    # --- Hydrogen ---
    h_s_valence = 3 + 1 + 1
    h_p_polarization = 1
    h_total = h_s_valence + h_p_polarization

    print("For a Hydrogen (H) atom:")
    # The prompt requires printing each number in the equation
    print(f"Primitives = (3+1+1) [s valence] + 1 [p polarization] = {h_s_valence} + {h_p_polarization} = {h_total}")
    print("-" * 40)

    # --- Second-Row Heavy Atom (e.g., Carbon) ---
    heavy_core = 6
    heavy_s_valence = 3 + 1 + 1
    heavy_p_valence = 3 + 1 + 1
    heavy_d_polarization = 1
    heavy_total = heavy_core + heavy_s_valence + heavy_p_valence + heavy_d_polarization

    print("For a second-row heavy atom (like Carbon, C):")
    # The prompt requires printing each number in the equation
    print(f"Primitives = 6 [core] + (3+1+1) [s valence] + (3+1+1) [p valence] + 1 [d polarization]")
    print(f"           = {heavy_core} + {heavy_s_valence} + {heavy_p_valence} + {heavy_d_polarization} = {heavy_total}")
    print("-" * 40)

    # --- Example Molecule: Methane (CH4) ---
    print("As a full example, let's calculate the total for Methane (CH4):\n")
    num_c = 1
    num_h = 4

    total_for_c = num_c * heavy_total
    total_for_h = num_h * h_total
    final_total = total_for_c + total_for_h

    print("Total Primitives = (Num C atoms * Primitives per C) + (Num H atoms * Primitives per H)")
    print(f"Total Primitives = ({num_c} * {heavy_total}) + ({num_h} * {h_total})")
    print(f"Total Primitives = {total_for_c} + {total_for_h}")
    print(f"Total Primitives = {final_total}")
    
    # Return the final calculated value for the example molecule
    return final_total

# Execute the function and capture the final answer for the example
final_answer = calculate_and_print_primitives()

# The final answer is the total number of primitive Gaussians for Methane (CH4).
# print(f"\n<<<For CH4: {final_answer}>>>") # This is for internal check, not final output format.