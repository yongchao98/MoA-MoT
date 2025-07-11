def find_synonym_indices():
    """
    This function stores and prints the indices of species from the 1872 survey
    that are considered synonyms in modern taxonomy.
    """
    # Indices of species identified as synonyms based on taxonomic research.
    synonym_indices = [4, 5, 6, 7, 8, 9, 11, 14, 15, 18, 19]

    # Sort the indices in ascending order (they are already sorted in this case).
    synonym_indices.sort()

    # Convert each integer index to a string.
    indices_as_strings = map(str, synonym_indices)

    # Join the string representations with a comma.
    final_output = ",".join(indices_as_strings)

    # Print the final formatted string.
    print(final_output)

find_synonym_indices()