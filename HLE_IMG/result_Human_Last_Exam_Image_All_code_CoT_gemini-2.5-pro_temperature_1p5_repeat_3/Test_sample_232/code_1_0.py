def identify_disorder():
    """
    Identifies the atrophied organ part and the associated disorder from the provided image.
    The image shows marked atrophy of the caudate nucleus, which is a classic sign of Huntington's disease.
    """
    # The atrophied part of the organ (brain)
    atrophied_structure = "caudate nucleus"

    # The name of the disorder
    disorder = "Huntington's disease"

    # Print the answer in the required "organ, disorder" format
    print(f"{atrophied_structure}, {disorder}")

identify_disorder()