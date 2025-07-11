def calculate_molecular_weight():
    """
    Calculates the molecular weight of the products A and B, which are isomers
    with the molecular formula C19H22O3.
    """
    # Atomic weights (g/mol)
    atomic_weight = {
        'C': 12.011,
        'H': 1.008,
        'O': 15.999
    }

    # Molecular formula for both products A and B
    formula = {'C': 19, 'H': 22, 'O': 3}

    # Calculate molecular weight
    mw_carbon = formula['C'] * atomic_weight['C']
    mw_hydrogen = formula['H'] * atomic_weight['H']
    mw_oxygen = formula['O'] * atomic_weight['O']
    
    total_mw = mw_carbon + mw_hydrogen + mw_oxygen

    print("The molecular formula for both products A and B is C19H22O3.")
    print("The calculation of the molecular weight is as follows:")
    print(f"({formula['C']} * {atomic_weight['C']}) + ({formula['H']} * {atomic_weight['H']}) + ({formula['O']} * {atomic_weight['O']})")
    print(f"= {mw_carbon:.3f} + {mw_hydrogen:.3f} + {mw_oxygen:.3f}")
    print(f"= {total_mw:.3f} g/mol")

if __name__ == "__main__":
    calculate_molecular_weight()
