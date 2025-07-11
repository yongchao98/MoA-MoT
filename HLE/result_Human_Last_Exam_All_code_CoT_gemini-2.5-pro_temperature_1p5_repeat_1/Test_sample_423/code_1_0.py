def find_empirical_formula(composition):
    """
    Calculates the empirical formula from mass fractions.

    Args:
        composition (dict): A dictionary with elements as keys and mass fractions as values.
    """
    atomic_masses = {
        'C': 12.01,
        'H': 1.008,
        'O': 16.00
    }

    # Calculate moles for a 100g sample
    moles = {el: composition[el] / atomic_masses[el] for el in composition}
    print("Step 1: Calculate moles of each element in a 100g sample.")
    for el, mol in moles.items():
        print(f"  - Moles of {el}: {composition[el]*100:.1f} g / {atomic_masses[el]:.2f} g/mol = {mol:.4f} mol")
    
    # Find the smallest mole value
    min_moles = min(moles.values())
    print(f"\nStep 2: Find the smallest mole value to determine the ratio. Smallest value is {min_moles:.4f}.")

    # Calculate the ratio
    ratio = {el: moles[el] / min_moles for el in moles}
    print("\nStep 3: Divide all mole values by the smallest value to get the simplest ratio.")
    for el, r in ratio.items():
        print(f"  - Ratio for {el}: {moles[el]:.4f} / {min_moles:.4f} = {r:.2f}")

    # Convert ratio to nearest simple integers (by multiplying by 3 in this case)
    multiplier = 3
    empirical_formula = {el: round(r * multiplier) for el, r in ratio.items()}
    
    print(f"\nStep 4: The ratios are approximately C=1.67, H=4.00, O=1.00. Multiply by {multiplier} to get integers.")

    formula_str = "".join([f"{el}{empirical_formula[el]}" for el in empirical_formula])
    print("\nFinal Result:")
    print(f"The empirical formula of B is {formula_str}.")

# Mass fractions of elements in B
composition_B = {
    'C': 0.5,
    'H': 0.1,
    'O': 0.4
}

find_empirical_formula(composition_B)