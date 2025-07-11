def solve():
    """
    This function identifies the indices of species that perform photochemical synthesis from a predefined list.
    """
    # The indices of the species that perform photochemical synthesis are:
    # 2: Gloeochaete wittrockiana (alga)
    # 5: Halapricum salinum (archaeon)
    # 6: Aphanothece castagnei (cyanobacterium)
    # 7: Baileya pleniradiata (plant)
    # 9: Ectothiorhodosinus mongolicus (purple sulfur bacterium)
    # 10: Chlorobaculum tepidum (green sulfur bacterium)
    # 12: Gemmatimonas phototrophica (phototrophic bacterium)
    
    correct_indices = [2, 5, 6, 7, 9, 10, 12]
    
    # Convert each integer index to a string
    indices_as_strings = [str(i) for i in correct_indices]
    
    # Join the string representations with a comma
    result = ",".join(indices_as_strings)
    
    print(result)

solve()