def solve_graphene_classification():
    """
    Analyzes and classifies three graphene nanoribbon band structures.

    The analysis is based on visual inspection of the provided image:
    - Edge Type (A=Armchair, Z=Zigzag): All three plots show band structures
      symmetric around k=0, characteristic of Armchair nanoribbons.
    - Width (N): Determined by counting the number of conduction bands (E>0).
    - Band Type (0=metallic, 1=semiconducting): Determined by the presence
      or absence of a band gap at the Fermi level (E=0).

    Analysis Results:
    - Ribbon 1: Edge='A', Width=8 (8 bands above E=0), Band=1 (semiconducting gap).
    - Ribbon 2: Edge='A', Width=5 (5 bands above E=0), Band=0 (metallic, no gap).
    - Ribbon 3: Edge='A', Width=6 (6 bands above E=0), Band=1 (semiconducting gap).
    """

    # Store the classification data for each ribbon
    # Format: (Edge_Type, Width, Band_Type)
    classifications = [
        ('A', 8, 1),  # Ribbon 1
        ('A', 5, 0),  # Ribbon 2
        ('A', 6, 1)   # Ribbon 3
    ]

    # Build the final classification string by concatenating each part
    final_string_parts = []
    for edge, width, band in classifications:
        # The prompt requires outputting each number.
        # We build the string from these numbers.
        part = f"{edge}{width}{band}"
        final_string_parts.append(part)

    final_classification_string = "".join(final_string_parts)

    # Print the final concatenated string
    print(final_classification_string)

solve_graphene_classification()