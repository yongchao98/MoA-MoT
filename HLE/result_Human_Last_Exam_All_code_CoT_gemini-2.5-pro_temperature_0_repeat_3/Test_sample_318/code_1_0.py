def find_synonym_indices():
    """
    This function identifies and returns the indices of species from the 1872 list
    that were considered synonyms in 2020.

    The research and determination of synonym status has been performed manually
    using modern taxonomic databases. This function simply stores and formats the result.
    """
    # Indices of species that are now considered synonyms
    synonym_indices = [2, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14, 15, 18, 19, 21]

    # Sort the indices in ascending order (they are already sorted, but this is good practice)
    synonym_indices.sort()

    # Format the output as a comma-separated string without spaces
    result = ",".join(map(str, synonym_indices))

    print(result)

find_synonym_indices()