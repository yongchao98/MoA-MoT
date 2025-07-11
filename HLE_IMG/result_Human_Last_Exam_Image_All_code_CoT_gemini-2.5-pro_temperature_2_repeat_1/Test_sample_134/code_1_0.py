def name_molecule():
    """
    This script identifies the components of the given molecule, calculates its
    molecular formula, and prints its name.
    """
    components = {
        'p-phenylene': {'formula': (6, 4), 'count': 6},
        'o-phenylene': {'formula': (6, 4), 'count': 2},
        'phenanthrylene': {'formula': (14, 8), 'count': 2},
        'ethynylene': {'formula': (2, 0), 'count': 8},
        'vinylene': {'formula': (2, 2), 'count': 2}
    }

    total_C = 0
    total_H = 0

    carbon_calculation_str = []
    hydrogen_calculation_str = []
    
    print("Calculating the molecular formula based on the molecule's components:\n")
    
    # Sort components for a consistent output order
    sorted_components = sorted(components.items())

    for name, data in sorted_components:
        C, H = data['formula']
        count = data['count']
        
        c_contrib = C * count
        h_contrib = H * count
        
        total_C += c_contrib
        total_H += h_contrib
        
        carbon_calculation_str.append(str(c_contrib))
        hydrogen_calculation_str.append(str(h_contrib))
        
        print(f"- {name.replace('_', ' ').capitalize()} (C{C}H{H}): {count} units")
        print(f"  Carbon: {count} * {C} = {c_contrib}")
        print(f"  Hydrogen: {count} * {H} = {h_contrib}")
        print("-" * 20)

    print("\nTotal Carbon calculation:")
    print(f"C = {' + '.join(carbon_calculation_str)} = {total_C}")
    
    print("\nTotal Hydrogen calculation:")
    print(f"H = {' + '.join(hydrogen_calculation_str)} = {total_H}")
    
    print("\n" + "="*30)
    print(f"Final Molecular Formula: C{total_C}H{total_H}")
    print("Molecule Name: Warped Carbon Nanophene C8 (WCNP-C8)")
    print("="*30)

name_molecule()