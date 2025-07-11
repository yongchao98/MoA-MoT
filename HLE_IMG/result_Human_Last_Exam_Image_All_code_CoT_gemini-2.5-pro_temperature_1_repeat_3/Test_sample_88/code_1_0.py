def calculate_molecular_weight():
    """
    Calculates the molecular weight of Product C (C11H16N2O3) and
    prints the detailed calculation.
    """
    # Standard atomic weights of the elements
    atomic_weights = {
        'C': 12.011,
        'H': 1.008,
        'N': 14.007,
        'O': 15.999
    }

    # Molecular formula of Product C
    formula_C = {'C': 11, 'H': 16, 'N': 2, 'O': 3}

    # Unpack the counts for clarity
    c_count = formula_C['C']
    h_count = formula_C['H']
    n_count = formula_C['N']
    o_count = formula_C['O']

    # Calculate the weight contribution of each element
    c_weight = c_count * atomic_weights['C']
    h_weight = h_count * atomic_weights['H']
    n_weight = n_count * atomic_weights['N']
    o_weight = o_count * atomic_weights['O']

    # Calculate the total molecular weight
    total_mw = c_weight + h_weight + n_weight + o_weight

    print("Molecular weight calculation for Product C (C11H16N2O3):")
    print("=======================================================")
    print(f"The formula for calculation is: ({c_count} * C) + ({h_count} * H) + ({n_count} * N) + ({o_count} * O)")
    print("\nBreaking it down:")
    
    # Print each part of the equation
    print(f"Contribution from Carbon (C):   {c_count} * {atomic_weights['C']} = {c_weight:.3f}")
    print(f"Contribution from Hydrogen (H): {h_count} * {atomic_weights['H']}  = {h_weight:.3f}")
    print(f"Contribution from Nitrogen (N):  {n_count} * {atomic_weights['N']}  = {n_weight:.3f}")
    print(f"Contribution from Oxygen (O):    {o_count} * {atomic_weights['O']}  = {o_weight:.3f}")
    
    print("-------------------------------------------------------")
    # Print the final sum
    print("Final equation:")
    print(f"{c_weight:.3f} + {h_weight:.3f} + {n_weight:.3f} + {o_weight:.3f} = {total_mw:.3f}")
    print("=======================================================")
    print(f"The molecular weight of Product C is {total_mw:.3f} g/mol.")


if __name__ == "__main__":
    calculate_molecular_weight()