def calculate_product_properties():
    """
    Calculates and prints the molecular formula and weight of the reaction product.

    The reaction is an anionic Oxy-Cope rearrangement, which is an isomerization.
    Therefore, the product has the same molecular formula as the starting material.
    """

    # --- Step 1: Define Atomic Weights ---
    atomic_weights = {
        'C': 12.011,  # Carbon
        'H': 1.008,   # Hydrogen
        'O': 15.999,  # Oxygen
        'Si': 28.085  # Silicon
    }

    # --- Step 2: Define Molecular Formula of the Product ---
    # The reaction is an isomerization, so the product formula is the same as the starting material's.
    # Starting Material: (1S,2R,4S)-2-((S)-4-((tert-butyldimethylsilyl)oxy)cyclopent-1-en-1-yl)-7,7-dimethoxybicyclo[2.2.1]hept-5-en-2-ol
    # Formula: C20 H34 O4 Si
    formula = {
        'C': 20,
        'H': 34,
        'O': 4,
        'Si': 1
    }

    # --- Step 3: Explain the Product and Calculate Molecular Weight ---
    print("The reaction is an anionic Oxy-Cope rearrangement, which is a type of isomerization.")
    print("The product is a fused bicyclic ketone with the same molecular formula as the starting material.")
    print("-" * 30)
    print(f"Product Molecular Formula: C{formula['C']}H{formula['H']}O{formula['O']}Si")
    print("-" * 30)
    
    # Calculate the molecular weight
    molecular_weight = (formula['C'] * atomic_weights['C'] +
                        formula['H'] * atomic_weights['H'] +
                        formula['O'] * atomic_weights['O'] +
                        formula['Si'] * atomic_weights['Si'])

    # --- Step 4: Print the detailed calculation ---
    print("Molecular Weight Calculation:")
    print(
        f"({formula['C']} * {atomic_weights['C']}) + "
        f"({formula['H']} * {atomic_weights['H']}) + "
        f"({formula['O']} * {atomic_weights['O']}) + "
        f"({formula['Si']} * {atomic_weights['Si']}) = {molecular_weight:.3f} g/mol"
    )

# Execute the function
calculate_product_properties()