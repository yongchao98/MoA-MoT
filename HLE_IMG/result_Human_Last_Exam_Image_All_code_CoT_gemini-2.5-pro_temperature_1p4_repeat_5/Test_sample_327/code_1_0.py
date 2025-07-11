def solve_graphene_puzzle():
    """
    This function analyzes the properties of three graphene nanoribbon band structures
    and generates a concatenated classification string.
    """

    # --- Analysis of each plot based on visual inspection ---

    # Plot 1 Analysis
    edge_1 = 'A'  # Armchair: No flat band at E=0, symmetric about k=0.
    width_1 = 8   # Counted 8 bands above the Fermi level (E=0).
    band_1 = 1    # Semiconducting: Clear band gap at E=0.
    
    # Plot 2 Analysis
    edge_2 = 'Z'  # Zigzag: Prominent flat band at E=0.
    width_2 = 7   # Counted 7 bands above the Fermi level.
    band_2 = 0    # Metallic: Bands touch at E=0, no gap.

    # Plot 3 Analysis
    edge_3 = 'A'  # Armchair: No flat band at E=0, symmetric about k=0.
    width_3 = 6   # Counted 6 bands above the Fermi level.
    band_3 = 1    # Semiconducting: Clear band gap at E=0.

    # --- Constructing the output string ---

    # Create the classification string for each plot
    class_1 = f"{edge_1}{width_1}{band_1}"
    class_2 = f"{edge_2}{width_2}{band_2}"
    class_3 = f"{edge_3}{width_3}{band_3}"

    # Concatenate the individual strings
    final_classification = class_1 + class_2 + class_3
    
    # Print the final result
    print("Classification for Plot 1:")
    print(f"Edge: {edge_1}, Width: {width_1}, Band: {band_1} -> {class_1}")
    print("\nClassification for Plot 2:")
    print(f"Edge: {edge_2}, Width: {width_2}, Band: {band_2} -> {class_2}")
    print("\nClassification for Plot 3:")
    print(f"Edge: {edge_3}, Width: {width_3}, Band: {band_3} -> {class_3}")
    
    print("\nFinal Concatenated String:")
    print(final_classification)

solve_graphene_puzzle()