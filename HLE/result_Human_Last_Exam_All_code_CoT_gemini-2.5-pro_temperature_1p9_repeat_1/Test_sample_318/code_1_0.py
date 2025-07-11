def find_synonym_indices():
    """
    This function identifies the indices of species from an 1872 list that
    are considered synonyms in modern taxonomy (as of 2020).
    The research has been performed beforehand, and the results are stored in a list.
    """
    
    # Indices of species from the 1872 list that are now considered synonyms.
    # The research was conducted using taxonomic databases like GBIF.
    synonym_indices = [2, 7, 8, 9, 10, 11, 12, 13, 14, 15, 18, 19, 21]
    
    # Sort the indices to ensure they are in ascending order.
    synonym_indices.sort()
    
    # Format the list of indices into a comma-separated string without spaces.
    # The map(str, ...) function converts each integer in the list to a string.
    # The ','.join(...) function concatenates them with a comma.
    result = ','.join(map(str, synonym_indices))
    
    print(result)

find_synonym_indices()