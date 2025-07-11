def calculate_molecular_weight():
    """
    Calculates the molecular weight of the starting material,
    (3,4-dihydro-2H-pyrrol-5-yl)proline, with the formula C9H14N2O2.
    """
    # Atomic weights of the elements
    atomic_weights = {
        'C': 12.011,
        'H': 1.008,
        'N': 14.007,
        'O': 15.999
    }

    # Molecular formula of the starting material
    formula = {
        'C': 9,
        'H': 14,
        'N': 2,
        'O': 2
    }

    # Calculate the molecular weight
    mw = (formula['C'] * atomic_weights['C'] +
          formula['H'] * atomic_weights['H'] +
          formula['N'] * atomic_weights['N'] +
          formula['O'] * atomic_weights['O'])

    # Print the equation as requested
    print(f"The molecular formula of the starting material is C9H14N2O2.")
    print(f"The molecular weight is calculated as follows:")
    print(
        f"{formula['C']} * {atomic_weights['C']} (C) + "
        f"{formula['H']} * {atomic_weights['H']} (H) + "
        f"{formula['N']} * {atomic_weights['N']} (N) + "
        f"{formula['O']} * {atomic_weights['O']} (O) = "
        f"{mw:.3f}"
    )
    
    # Returning the final value for the answer block
    return mw

if __name__ == '__main__':
    molecular_weight = calculate_molecular_weight()
    # The final answer in the requested format
    # print(f"<<<{molecular_weight:.3f}>>>")