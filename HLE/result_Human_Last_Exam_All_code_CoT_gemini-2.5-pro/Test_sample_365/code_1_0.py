def calculate_molecular_formula():
    """
    Calculates and prints the molecular formula for the product of the reaction.
    The reaction is an intramolecular rearrangement (anionic oxy-Cope), so the product
    is an isomer of the starting material, having the same molecular formula.
    """
    # Atom counts for the starting material:
    # (1S,2R,4S)-2-((S)-4-((tert-butyldimethylsilyl)oxy)cyclopent-1-en-1-yl)-7,7-dimethoxybicyclo[2.2.1]hept-5-en-2-ol
    # Norbornene core (C7) + dimethoxy (C2) = C9
    # Cyclopentenyl substituent (C5) = C5
    # TBDMS group (Si(CH3)2(C(CH3)3)) = C6
    # Total C = 7 + 2 + 5 + 6 = 20
    atom_counts = {
        'C': 20,
        'H': 34,
        'O': 4,
        'Si': 1
    }

    print("The reaction is an anionic oxy-Cope rearrangement, which is an isomerization.")
    print("Therefore, the product has the same molecular formula as the starting material.")
    print("\nCalculating the molecular formula based on atom counts:")
    
    formula = ""
    # Standard order is C, H, then alphabetical for the rest.
    elements_order = ['C', 'H', 'O', 'Si']
    for element in elements_order:
        count = atom_counts.get(element, 0)
        if count > 0:
            print(f"Number of {element} atoms: {count}")
            formula += element
            if count > 1:
                formula += str(count)

    print("\n---")
    print(f"The resulting molecular formula of the product is: {formula}")
    print("---")
    print("\nThe product is a bicyclic ketone, formed by converting the tertiary alcohol to a ketone")
    print("and expanding the original bicyclo[2.2.1]heptene ring system.")
    print("The silyl ether and ketal protecting groups remain intact.")

calculate_molecular_formula()