def calculate_bragg_reflections():
    """
    Calculates and prints the number of Bragg reflections for {200}, {220}, and {222}
    families of planes in a rhombohedrally distorted perovskite (R3m space group).
    """

    # Due to the symmetry reduction from cubic to rhombohedral, the Bragg peaks split.
    # The number of resulting reflections for each family is determined by how the planes
    # in the family are partitioned into non-equivalent sets in the lower symmetry.

    # Number of reflections for the {200} family
    reflections_200 = 2

    # Number of reflections for the {220} family
    reflections_220 = 2

    # Number of reflections for the {222} family
    reflections_222 = 2

    # Calculate the total number of reflections from these families
    total_reflections = reflections_200 + reflections_220 + reflections_222

    print("In a rhombohedral (R3m) crystal system, the pseudocubic reflections split as follows:")
    print("-" * 75)
    print(f"The {200} family of planes splits into {reflections_200} Bragg reflections.")
    print(f"The {220} family of planes splits into {reflections_220} Bragg reflections.")
    print(f"The {222} family of planes splits into {reflections_222} Bragg reflections.")
    print("-" * 75)
    
    # The final "equation" showing the breakdown
    print("\nThe equation representing the total count of these observed reflections is:")
    print(f"{reflections_200} (from {{200}}) + {reflections_220} (from {{220}}) + {reflections_222} (from {{222}}) = {total_reflections} reflections")

# Execute the function to print the results
calculate_bragg_reflections()