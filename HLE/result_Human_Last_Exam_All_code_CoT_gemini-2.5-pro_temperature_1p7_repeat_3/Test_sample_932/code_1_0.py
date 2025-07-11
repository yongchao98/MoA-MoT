def find_unlikely_tribes():
    """
    Identifies the indices of insect tribes whose immatures are unlikely to be
    collected using a beat-sheet method.

    The analysis of each tribe is as follows:
    1) Apis: Immatures live inside a hive. Unlikely.
    2) Melipotini: Immatures (caterpillars) live on foliage. Likely.
    3) Eupholini: Immatures (larvae) are internal plant borers. Unlikely.
    4) Acritini: Immatures live in dung, carrion, or soil. Unlikely.
    5) Oxyptilini: Immatures (caterpillars) live on plants. Likely.
    6) Dictyophorini: Immatures (nymphs) live on foliage. Likely.
    7) Acanthocerini: Immatures (grubs) live in soil or decaying wood. Unlikely.

    The function will print the indices of the unlikely tribes.
    """
    unlikely_indices = [1, 3, 4, 7]
    
    # Sort the indices in ascending order (they are already sorted but this is good practice)
    unlikely_indices.sort()
    
    # Convert numbers to strings to join them with commas
    # The final format should be "1,3,4,7"
    output_string = ",".join(map(str, unlikely_indices))
    
    print(output_string)

find_unlikely_tribes()