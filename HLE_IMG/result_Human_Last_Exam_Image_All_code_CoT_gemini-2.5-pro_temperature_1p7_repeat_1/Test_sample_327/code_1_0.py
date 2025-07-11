def solve_graphene_classification():
    """
    This function classifies three graphene nanoribbon band structures
    based on their edge type, width, and electronic properties, then
    concatenates the results into a single string.
    """

    # --- Analysis of Ribbon 1 ---
    # Edge type: The band structure shows a symmetric gap at k=0, which is characteristic of an Armchair (A) nanoribbon.
    edge_1 = 'A'
    # Width: Counting the bands above the E=0 axis (conduction bands) gives 8. The problem format suggests the width parameter
    # is the number of atoms in the unit cell, which is 2 * (number of dimer lines).
    num_dimer_lines_1 = 8
    width_1 = 2 * num_dimer_lines_1
    # Band type: There is a clear energy gap around the Fermi level (E=0), so it is semiconducting (1).
    band_1 = 1
    class_1 = f"{edge_1}{width_1}{band_1}"

    print("Classification for Ribbon 1:")
    print(f"Edge='{edge_1}', N_dimer_lines={num_dimer_lines_1}, Width_param={width_1}, Band={band_1} => {class_1}")

    # --- Analysis of Ribbon 2 ---
    # Edge type: The presence of flat bands near E=0 is a signature of a Zigzag (Z) nanoribbon.
    edge_2 = 'Z'
    # Width: Counting the conduction bands gives 7. The width parameter is 2 * (number of zigzag chains).
    num_zigzag_chains_2 = 7
    width_2 = 2 * num_zigzag_chains_2
    # Band type: The bands cross the Fermi level, indicating it is metallic (0).
    band_2 = 0
    class_2 = f"{edge_2}{width_2}{band_2}"
    
    print("\nClassification for Ribbon 2:")
    print(f"Edge='{edge_2}', N_zigzag_chains={num_zigzag_chains_2}, Width_param={width_2}, Band={band_2} => {class_2}")

    # --- Analysis of Ribbon 3 ---
    # Edge type: Similar to Ribbon 1, the gap at k=0 indicates an Armchair (A) nanoribbon.
    edge_3 = 'A'
    # Width: Counting the conduction bands gives 9.
    num_dimer_lines_3 = 9
    width_3 = 2 * num_dimer_lines_3
    # Band type: It is semiconducting (1) due to the visible band gap.
    band_3 = 1
    class_3 = f"{edge_3}{width_3}{band_3}"

    print("\nClassification for Ribbon 3:")
    print(f"Edge='{edge_3}', N_dimer_lines={num_dimer_lines_3}, Width_param={width_3}, Band={band_3} => {class_3}")

    # --- Final Concatenated String ---
    final_classification = class_1 + class_2 + class_3
    print("\nFinal concatenated string:")
    print(final_classification)

solve_graphene_classification()