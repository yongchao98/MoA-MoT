def calculate_molecular_weight():
    """
    Identifies Compound 1 and calculates its molecular formula and weight.
    
    The reaction described is the formation of an O-allyl thionocarbonate,
    which immediately undergoes a [3,3]-sigmatropic rearrangement (thiono-Claisen
    rearrangement) to form a more stable S-allyl thiocarbonate. This rearrangement
    explains the change in the NMR spectrum.

    Compound 1 is: O-(p-tolyl) S-(3,7-dimethylocta-1,6-dien-3-yl) thiocarbonate.

    This script calculates its molecular weight.
    """

    # Molecular formula of Compound 1 is C18H24O2S
    formula = {'C': 18, 'H': 24, 'O': 2, 'S': 1}
    
    # Standard atomic weights
    atomic_weights = {
        'C': 12.011,
        'H': 1.008,
        'O': 15.999,
        'S': 32.06
    }
    
    print("Compound 1 is the rearranged product: O-(p-tolyl) S-(3,7-dimethylocta-1,6-dien-3-yl) thiocarbonate.")
    print("Its molecular formula is C18H24O2S.")
    print("\nCalculating the molecular weight:\n")

    total_weight = 0
    calculation_steps = []

    # Iterating through the formula to build the output string and calculate the total weight
    for element, count in formula.items():
        weight = atomic_weights[element]
        element_total_weight = count * weight
        total_weight += element_total_weight
        # "output each number in the final equation!"
        print(f"Contribution from {element}: {count} * {weight} = {element_total_weight:.3f}")
        calculation_steps.append(f"({count} * {weight})")

    print("\nFinal Molecular Weight Calculation:")
    equation_str = " + ".join(calculation_steps)
    print(f"{equation_str} = {total_weight:.3f} g/mol")
    
    print(f"\nThe molecular weight of Compound 1 is {total_weight:.3f} g/mol.")

# Execute the function
calculate_molecular_weight()