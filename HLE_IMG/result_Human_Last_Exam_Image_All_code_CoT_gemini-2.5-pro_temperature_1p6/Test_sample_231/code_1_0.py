def calculate_molecular_weight():
    """
    Calculates and prints the molecular weight of compound C (C23H23NO8).
    The derivation of the molecular formula is based on the reaction scheme provided.
    The function prints the contribution of each element to the total molecular weight.
    """
    # Atomic weights of the constituent elements (g/mol)
    atomic_weights = {
        'C': 12.011,
        'H': 1.008,
        'N': 14.007,
        'O': 15.999
    }

    # Molecular formula of the final neutral compound C
    molecular_formula = {
        'C': 23,
        'H': 23,
        'N': 1,
        'O': 8
    }

    total_mw = 0.0

    print("Step-by-step calculation of the molecular weight of Compound C (C23H23NO8):")
    print("-" * 65)

    # Carbon contribution
    num_C = molecular_formula['C']
    aw_C = atomic_weights['C']
    mass_C = num_C * aw_C
    total_mw += mass_C
    print(f"Contribution from Carbon (C):   {num_C} atoms * {aw_C} u = {mass_C:.3f} u")

    # Hydrogen contribution
    num_H = molecular_formula['H']
    aw_H = atomic_weights['H']
    mass_H = num_H * aw_H
    total_mw += mass_H
    print(f"Contribution from Hydrogen (H): {num_H} atoms * {aw_H} u = {mass_H:.3f} u")

    # Nitrogen contribution
    num_N = molecular_formula['N']
    aw_N = atomic_weights['N']
    mass_N = num_N * aw_N
    total_mw += mass_N
    print(f"Contribution from Nitrogen (N): {num_N} atom  * {aw_N} u = {mass_N:.3f} u")

    # Oxygen contribution
    num_O = molecular_formula['O']
    aw_O = atomic_weights['O']
    mass_O = num_O * aw_O
    total_mw += mass_O
    print(f"Contribution from Oxygen (O):   {num_O} atoms * {aw_O} u = {mass_O:.3f} u")

    print("-" * 65)
    print("Final equation for the total molecular weight:")
    print(f"Total MW = {mass_C:.3f} (from C) + {mass_H:.3f} (from H) + {mass_N:.3f} (from N) + {mass_O:.3f} (from O)")
    print(f"Total Molecular Weight of Compound C = {total_mw:.3f} u")

if __name__ == '__main__':
    calculate_molecular_weight()