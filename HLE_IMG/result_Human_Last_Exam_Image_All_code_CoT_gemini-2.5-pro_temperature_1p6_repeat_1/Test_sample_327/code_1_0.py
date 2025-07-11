def classify_nanoribbons():
    """
    This function classifies three graphene nanoribbon band structures based on visual analysis
    and prints the concatenated result.

    Classification Rules:
    - Edge: 'A' for Armchair, 'Z' for Zigzag.
    - Width (N): Number of dimer lines (equals the number of conduction bands).
    - Band: '1' for semiconducting (has a band gap), '0' for metallic (no band gap).
    """

    # --- Classification based on visual analysis of the plots ---

    # Ribbon 1 Classification
    # Edge: Armchair (symmetric, no flat edge state at E=0)
    # Width: 8 bands counted above E=0
    # Band: Semiconducting (clear band gap at E=0)
    r1_edge = 'A'
    r1_width = 8
    r1_band = 1

    # Ribbon 2 Classification
    # Edge: Armchair (symmetric, no flat edge state at E=0)
    # Width: 5 bands counted above E=0
    # Band: Metallic (bands touch at E=0)
    r2_edge = 'A'
    r2_width = 5
    r2_band = 0

    # Ribbon 3 Classification
    # Edge: Zigzag (characteristic flat edge state at E=0)
    # Width: 6 bands counted above E=0
    # Band: Metallic (no band gap)
    r3_edge = 'Z'
    r3_width = 6
    r3_band = 0

    # --- Constructing and printing the final output ---

    # Individual classifications
    class1 = f"{r1_edge}{r1_width}{r1_band}"
    class2 = f"{r2_edge}{r2_width}{r2_band}"
    class3 = f"{r3_edge}{r3_width}{r3_band}"

    # Final concatenated string
    final_result = class1 + class2 + class3
    
    print("The final classification string is built by concatenating the properties of each ribbon.")
    print(f"Ribbon 1: Edge='{r1_edge}', Width={r1_width}, Band={r1_band} -> {class1}")
    print(f"Ribbon 2: Edge='{r2_edge}', Width={r2_width}, Band={r2_band} -> {class2}")
    print(f"Ribbon 3: Edge='{r3_edge}', Width={r3_width}, Band={r3_band} -> {class3}")
    print("\nConcatenated classification:")
    print(final_result)

# Run the classification
classify_nanoribbons()