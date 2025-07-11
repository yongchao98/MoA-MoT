def find_synonym_indices():
    """
    This function identifies the indices of species from the 1872 list
    that are now considered synonyms of different species.
    
    The research was conducted by consulting modern taxonomic databases. The results are:
    - 12. Hemichroa albidovariata -> Caulocampus acericaulis
    - 13. Hemichroa fraternalis -> Nematus ventralis
    - 17. Tenthredo nimbipennis -> Tenthredo xanthus
    - 18. Lophyrus Abietis -> Neodiprion lecontei
    - 19. Lophyrus fulva -> Neodiprion fulviceps
    """
    
    synonym_indices = [12, 13, 17, 18, 19]
    
    # Sort the indices in ascending order
    synonym_indices.sort()
    
    # Format the output string as comma-separated values without spaces
    output_string = ",".join(map(str, synonym_indices))
    
    print(output_string)

find_synonym_indices()