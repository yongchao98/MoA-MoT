def calculate_primitives_6_311g_star_star(symbol: str):
    """
    Calculates the number of primitive Gaussian functions in the 6-311G**
    basis set for a given element (H to Ar).

    Args:
        symbol: The atomic symbol of the element (e.g., 'H', 'C', 'Si').
    """
    s_valence_prims = 3 + 1 + 1
    p_valence_prims = 3 + 1 + 1

    # Atomic numbers for elements from H to Ar
    atomic_numbers = {
        'H': 1, 'HE': 2, 'LI': 3, 'BE': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'NE': 10,
        'NA': 11, 'MG': 12, 'AL': 13, 'SI': 14, 'P': 15, 'S': 16, 'CL': 17, 'AR': 18
    }

    z = atomic_numbers.get(symbol.upper())

    if z is None:
        print(f"Sorry, the element '{symbol}' is not supported by this script (supported: H-Ar).")
        return

    print(f"Calculating primitive Gaussians for {symbol} (Z={z}) in 6-311G** basis set:")
    print("-" * 20)

    # Period 1: H, He
    if 1 <= z <= 2:
        # For H and He, 1s is the valence orbital.
        # The second '*' adds p-polarization functions.
        pol_prims = 1
        total = s_valence_prims + pol_prims
        print("Valence 1s contribution: 3 + 1 + 1 = 5")
        print("Polarization 'p' contribution: 1")
        print(f"\nFinal Equation: (3+1+1) + 1 = {total}")

    # Period 2: Li - Ne
    elif 3 <= z <= 10:
        # Core is 1s orbital, described by 6 primitives.
        core_prims = 6
        # Valence are 2s and 2p orbitals.
        # The first '*' adds d-polarization functions.
        pol_prims = 1
        total = core_prims + s_valence_prims + p_valence_prims + pol_prims
        print(f"Core 1s contribution: {core_prims}")
        print(f"Valence 2s contribution: 3 + 1 + 1 = {s_valence_prims}")
        print(f"Valence 2p contribution: 3 + 1 + 1 = {p_valence_prims}")
        print(f"Polarization 'd' contribution: {pol_prims}")
        print(f"\nFinal Equation: {core_prims} + ({3}+{1}+{1}) + ({3}+{1}+{1}) + {pol_prims} = {total}")

    # Period 3: Na - Ar
    elif 11 <= z <= 18:
        # Core are 1s, 2s, 2p orbitals, each described by 6 primitives.
        core_1s = 6
        core_2s = 6
        core_2p = 6
        core_prims_total = core_1s + core_2s + core_2p
        # Valence are 3s and 3p orbitals.
        # The first '*' adds d-polarization functions.
        pol_prims = 1
        total = core_prims_total + s_valence_prims + p_valence_prims + pol_prims
        print(f"Core 1s, 2s, 2p contribution: {core_1s} + {core_2s} + {core_2p} = {core_prims_total}")
        print(f"Valence 3s contribution: 3 + 1 + 1 = {s_valence_prims}")
        print(f"Valence 3p contribution: 3 + 1 + 1 = {p_valence_prims}")
        print(f"Polarization 'd' contribution: {pol_prims}")
        print(f"\nFinal Equation: ({core_1s}+{core_2s}+{core_2p}) + ({3}+{1}+{1}) + ({3}+{1}+{1}) + {pol_prims} = {total}")

if __name__ == '__main__':
    # You can change the element symbol here to get the count for another atom.
    element_symbol = input("Enter an element symbol (e.g., H, C, O, Si): ")
    calculate_primitives_6_311g_star_star(element_symbol)
