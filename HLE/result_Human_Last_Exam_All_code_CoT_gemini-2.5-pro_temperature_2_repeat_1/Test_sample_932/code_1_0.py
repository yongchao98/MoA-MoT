def solve():
    """
    This function identifies the indices of insect tribes whose immatures are
    unlikely to be collected by the beat-sheet method based on their life habits.
    """
    
    # Analysis:
    # 1) Apis: Immatures are in a hive. Unlikely.
    # 2) Melipotini: Caterpillars on leaves. Likely.
    # 3) Eupholini: Larvae are internal borers. Unlikely.
    # 4) Acritini: Larvae in decaying matter. Unlikely.
    # 5) Oxyptilini: Caterpillars on leaves. Likely.
    # 6) Dictyophorini: Nymphs on plants. Likely.
    # 7) Acanthocerini: Larvae in rotting wood. Unlikely.
    
    # Indices of the tribes unlikely to be collected.
    unlikely_indices = [1, 3, 4, 7]
    
    # Sort the indices in ascending order (they are already sorted).
    unlikely_indices.sort()
    
    # Format the output as a string with indices separated by commas.
    # The str() function converts each integer to a string.
    # The ",".join() method concatenates the list of strings with a comma.
    result = ",".join(str(index) for index in unlikely_indices)
    
    print(result)

solve()