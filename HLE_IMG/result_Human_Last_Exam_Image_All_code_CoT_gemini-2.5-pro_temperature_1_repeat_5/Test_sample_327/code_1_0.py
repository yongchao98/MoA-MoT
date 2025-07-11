def solve_nanoribbon_classification():
    """
    This function classifies three graphene nanoribbon band structures and
    prints the concatenated classification string.
    """

    # --- Classification for Ribbon 1 ---
    # Edge Type: Armchair (A) - Symmetric band structure with a gap at k=0.
    # Width (N): 8 - Number of bands counted in the conduction band (E > 0).
    # Band Type: Semiconducting (1) - A clear band gap exists at the Fermi level (E=0).
    edge_1 = "A"
    width_1 = 8
    band_1 = 1
    class_1 = f"{edge_1}{width_1}{band_1}"

    # --- Classification for Ribbon 2 ---
    # Edge Type: Armchair (A) - Symmetric band structure.
    # Width (N): 5 - Number of bands counted in the conduction band.
    # Band Type: Metallic (0) - Valence and conduction bands touch at the Fermi level.
    edge_2 = "A"
    width_2 = 5
    band_2 = 0
    class_2 = f"{edge_2}{width_2}{band_2}"

    # --- Classification for Ribbon 3 ---
    # Edge Type: Zigzag (Z) - Presence of a characteristic flat band at the Fermi level.
    # Width (N): 6 - Number of bands counted in the conduction band at k=0.
    # Band Type: Metallic (0) - The flat band at the Fermi level ensures it's metallic.
    edge_3 = "Z"
    width_3 = 6
    band_3 = 0
    class_3 = f"{edge_3}{width_3}{band_3}"

    # Concatenate all classifications into a single string
    final_result = class_1 + class_2 + class_3
    
    # The final equation is the concatenation of the individual classifications.
    # We print the final string which contains all the identified numbers.
    print(final_result)

solve_nanoribbon_classification()