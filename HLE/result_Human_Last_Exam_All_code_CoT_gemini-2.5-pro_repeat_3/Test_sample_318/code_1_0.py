def find_synonym_indices():
    """
    This function prints the indices of species from the 1872 list that are
    now considered synonyms, based on modern taxonomic research.
    The list of indices is predetermined from the research.
    """
    # The indices of species that have been identified as synonyms.
    synonym_indices = [2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 18, 19, 21]

    # The problem states "output each number in the final equation!".
    # As there is no equation, we will just print the final list of indices.
    # The list is already sorted, but we sort it to be robust.
    synonym_indices.sort()

    # Convert each integer index to a string.
    string_indices = [str(i) for i in synonym_indices]

    # Join the string indices with a comma, without spaces.
    final_output = ",".join(string_indices)

    print(final_output)

find_synonym_indices()