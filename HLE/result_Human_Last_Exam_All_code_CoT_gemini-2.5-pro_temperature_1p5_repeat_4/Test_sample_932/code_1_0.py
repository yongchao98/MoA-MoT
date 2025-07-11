def find_unlikely_tribes():
    """
    Identifies the indices of insect tribes whose immatures are unlikely to be
    collected using a beat-sheet method.

    The analysis is as follows:
    - 1) Apis: Larvae are in hives. Unlikely.
    - 2) Melipotini: Caterpillars on foliage. Likely.
    - 3) Eupholini: Wood-boring larvae. Unlikely.
    - 4) Acritini: Larvae in dung/carrion. Unlikely.
    - 5) Oxyptilini: Caterpillars on foliage. Likely.
    - 6) Dictyophorini: Nymphs on plants. Likely.
    - 7) Acanthocerini: Larvae in ant/termite nests. Unlikely.

    This function will print the indices of the unlikely tribes.
    """
    unlikely_indices = [1, 3, 4, 7]
    
    # Sort the indices in ascending order
    unlikely_indices.sort()
    
    # Format the output string as requested
    output_string = ", ".join(map(str, unlikely_indices))
    
    print(output_string)

find_unlikely_tribes()