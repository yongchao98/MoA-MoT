import sys

def calculate_primitives_6_311G_star_star():
    """
    Calculates the number of primitive Gaussian functions in a 6-311G** basis set
    for a given atom (H to Ar).
    """
    # Electron configurations: (core s-orbitals, core p-orbital sets, valence s-orbitals, valence p-orbital sets)
    # A "p-orbital set" refers to the group of (px, py, pz) orbitals.
    atom_configs = {
        # Row 1
        'H' : (0, 0, 1, 0), 'HE': (0, 0, 1, 0),
        # Row 2
        'LI': (1, 0, 1, 0), 'BE': (1, 0, 1, 0),
        'B' : (1, 0, 1, 1), 'C' : (1, 0, 1, 1), 'N' : (1, 0, 1, 1),
        'O' : (1, 0, 1, 1), 'F' : (1, 0, 1, 1), 'NE': (1, 0, 1, 1),
        # Row 3
        'NA': (2, 1, 1, 0), 'MG': (2, 1, 1, 0),
        'AL': (2, 1, 1, 1), 'SI': (2, 1, 1, 1), 'P' : (2, 1, 1, 1),
        'S' : (2, 1, 1, 1), 'CL': (2, 1, 1, 1), 'AR': (2, 1, 1, 1),
    }

    # Get user input
    try:
        symbol = input("Enter an atomic symbol (e.g., C for Carbon, H for Hydrogen): ").strip().upper()
        if symbol not in atom_configs:
            print(f"Error: Atom '{symbol}' is not supported or is an invalid symbol. Please choose from H to Ar.")
            sys.exit(1)
    except (EOFError, KeyboardInterrupt):
        print("\nCalculation cancelled.")
        sys.exit(0)


    config = atom_configs[symbol]
    num_core_s, num_core_p_sets, num_val_s, num_val_p_sets = config

    # --- Calculate primitives based on 6-311G** rules ---

    # Core orbitals: 6 primitives per orbital
    core_s_primitives = num_core_s * 6
    # Each p-set has 3 orbitals (px, py, pz)
    core_p_primitives = num_core_p_sets * 3 * 6

    # Valence orbitals: 3+1+1 = 5 primitives per orbital
    valence_s_primitives = num_val_s * 5
    valence_p_primitives = num_val_p_sets * 3 * 5

    # Polarization functions (**)
    # H, He get one p-set (3 primitives)
    if symbol in ['H', 'HE']:
        polarization_primitives = 3
        pol_type = "p-type"
    # Heavy atoms get one d-set (6 primitives for Cartesian d)
    else:
        polarization_primitives = 6
        pol_type = "d-type"

    total_primitives = (core_s_primitives + core_p_primitives +
                        valence_s_primitives + valence_p_primitives +
                        polarization_primitives)

    # --- Print the detailed breakdown ---
    print(f"\n--- Calculation for {symbol} with 6-311G** ---")
    
    equation_parts = []
    if core_s_primitives > 0:
        print(f"Core s-orbitals: {num_core_s} orbital(s) * 6 primitives = {core_s_primitives}")
        equation_parts.append(str(core_s_primitives))
    if core_p_primitives > 0:
        print(f"Core p-orbitals: {num_core_p_sets*3} orbital(s) * 6 primitives = {core_p_primitives}")
        equation_parts.append(str(core_p_primitives))
    if valence_s_primitives > 0:
        print(f"Valence s-orbitals: {num_val_s} orbital(s) * (3+1+1) primitives = {valence_s_primitives}")
        equation_parts.append(str(valence_s_primitives))
    if valence_p_primitives > 0:
        print(f"Valence p-orbitals: {num_val_p_sets*3} orbital(s) * (3+1+1) primitives = {valence_p_primitives}")
        equation_parts.append(str(valence_p_primitives))

    print(f"Polarization ({pol_type}): 1 set = {polarization_primitives} primitives")
    equation_parts.append(str(polarization_primitives))

    equation_str = " + ".join(equation_parts)
    print("\nFinal Equation:")
    print(f"Total Primitives = {equation_str} = {total_primitives}")


if __name__ == '__main__':
    calculate_primitives_6_311G_star_star()