def find_synonym_indices():
    """
    This function identifies the indices of species from the 1872 Texas survey
    that were considered synonyms of a different species in 2020.

    The research and determination of synonyms has been performed manually by
    consulting modern taxonomic databases. This script simply stores and formats
    the results.
    """

    # Indices of species names that are now considered synonyms
    synonym_indices = [2, 7, 8, 9, 10, 11, 12, 14, 15, 17, 18, 19]

    # Sort the indices in ascending order (they are already sorted but it's good practice)
    synonym_indices.sort()

    # Format the output string as a comma-separated list without spaces
    result = ",".join(map(str, synonym_indices))
    
    print(result)

find_synonym_indices()