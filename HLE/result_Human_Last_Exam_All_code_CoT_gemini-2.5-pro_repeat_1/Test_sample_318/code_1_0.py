def find_synonym_indices():
    """
    This function identifies and prints the indices of species from the 1872 survey
    that are considered synonyms in modern taxonomy.
    The indices have been pre-determined through taxonomic research.
    """
    # List of indices for species that are now considered synonyms.
    synonym_indices = [2, 4, 5, 6, 7, 8, 9, 10, 11, 14, 15, 18, 19]

    # The list is already in ascending order.
    # Convert each number to a string and join them with a comma.
    # This also satisfies the requirement to "output each number in the final equation".
    output_string = ",".join(map(str, synonym_indices))

    print(output_string)

find_synonym_indices()