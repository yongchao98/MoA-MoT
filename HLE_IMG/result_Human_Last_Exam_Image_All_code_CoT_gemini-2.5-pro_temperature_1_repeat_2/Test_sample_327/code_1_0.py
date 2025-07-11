def solve_nanoribbon_classification():
    """
    This function classifies three graphene nanoribbon band structures based on visual analysis
    and concatenates the results into a single string.
    """

    # --- Classification based on visual analysis of the plots ---

    # Plot 1 Classification
    edge_1 = 'A'  # Armchair: Parabolic dispersion around k=0
    width_1 = 8   # Width: 8 bands counted in the conduction band (E > 0)
    band_1 = 1    # Band: Semiconducting, as there is a clear band gap at E=0

    # Plot 2 Classification
    edge_2 = 'Z'  # Zigzag: Characteristic flat bands near the Fermi level
    width_2 = 7   # Width: 7 bands counted in the conduction band
    band_2 = 0    # Band: Metallic, as bands cross the Fermi level (no gap)

    # Plot 3 Classification
    edge_3 = 'A'  # Armchair: Parabolic dispersion around k=0
    width_3 = 5   # Width: 5 bands counted in the conduction band
    band_3 = 0    # Band: Metallic, as bands touch at E=0, k=0 (zero gap)

    # --- Construct the final string ---

    # Concatenate the classifications for each ribbon
    class_1 = f"{edge_1}{width_1}{band_1}"
    class_2 = f"{edge_2}{width_2}{band_2}"
    class_3 = f"{edge_3}{width_3}{band_3}"

    # Combine all classifications into the final result string
    final_result = f"{class_1}{class_2}{class_3}"

    # The problem statement mentions "output each number in the final equation".
    # This is interpreted as constructing the string from its numeric and character components.
    # The final_result string already does this implicitly.
    # For clarity, we can print the components.
    
    print(f"Ribbon 1: Edge={edge_1}, Width={width_1}, Band={band_1} -> {class_1}")
    print(f"Ribbon 2: Edge={edge_2}, Width={width_2}, Band={band_2} -> {class_2}")
    print(f"Ribbon 3: Edge={edge_3}, Width={width_3}, Band={band_3} -> {class_3}")
    print("\nFinal concatenated classification:")
    print(final_result)


solve_nanoribbon_classification()