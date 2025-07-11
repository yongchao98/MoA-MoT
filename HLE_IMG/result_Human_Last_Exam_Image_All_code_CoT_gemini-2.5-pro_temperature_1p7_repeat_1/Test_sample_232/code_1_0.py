def solve_medical_image_query():
    """
    This script identifies the atrophied structure and the disorder from the provided image
    and prints the answer in the specified format.
    """

    # The atrophied portion of the organ, based on visual analysis.
    # The image shows severe atrophy of the caudate nucleus.
    atrophied_structure = "caudate nucleus"

    # The disorder associated with this specific type of atrophy.
    # Atrophy of the caudate nucleus is a hallmark of Huntington's disease.
    disorder_name = "Huntington's disease"

    # Print the answer in the "structure, disorder" format.
    print(f"{atrophied_structure}, {disorder_name}")

solve_medical_image_query()