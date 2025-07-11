def solve_medical_image_query():
    """
    This function identifies the atrophied structure and the disorder shown in the provided medical image.
    The image shows a coronal section of a brain with marked atrophy of the caudate nucleus,
    a characteristic finding in Huntington's disease.
    """
    # The structure that has atrophied, which is a part of the brain.
    atrophied_structure = "caudate nucleus"

    # The disorder associated with this specific type of atrophy.
    disorder = "Huntington's disease"

    # Print the answer in the "organ, disorder" format.
    # Note: "caudate nucleus" is the specific atrophied part of the organ (brain).
    print(f"{atrophied_structure}, {disorder}")

solve_medical_image_query()