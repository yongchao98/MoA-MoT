def find_unlikely_tribes():
    """
    This function identifies the indices of insect tribes whose immatures
    are unlikely to be collected using a beat-sheet method.

    The beat-sheet method is effective for insects on external foliage.

    Evaluation:
    1) Apis: Immatures in hives. Unlikely.
    2) Melipotini: Caterpillars on foliage. Likely.
    3) Eupholini: Larvae are wood-borers. Unlikely.
    4) Acritini: Larvae in dung, carrion, or nests. Unlikely.
    5) Oxyptilini: Caterpillars on plants. Likely.
    6) Dictyophorini: Nymphs on plants. Likely.
    7) Acanthocerini: Larvae in nests or rotting wood. Unlikely.

    The indices of the unlikely tribes are collected and sorted.
    """
    unlikely_indices = [1, 3, 4, 7]
    
    # Sort the indices in ascending order (they are already sorted)
    unlikely_indices.sort()
    
    # Print the final result as a comma-separated string
    output = ", ".join(map(str, unlikely_indices))
    print(output)

find_unlikely_tribes()