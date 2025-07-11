def solve_task():
    """
    This function identifies and prints the indices of taxa that perform
    photochemical synthesis from a predefined list.
    """
    # The indices of the taxa that undergo photochemical synthesis are:
    # 2: Gloeochaete wittrockiana (Alga, oxygenic photosynthesis)
    # 5: Halapricum salinum (Archaea, retinal-based phototrophy)
    # 6: Aphanothece castagnei (Cyanobacterium, oxygenic photosynthesis)
    # 7: Baileya pleniradiata (Plant, oxygenic photosynthesis)
    # 9: Ectothiorhodosinus mongolicus (Bacterium, anoxygenic photosynthesis)
    # 10: Chlorobaculum tepidum (Bacterium, anoxygenic photosynthesis)
    # 12: Gemmatimonas phototrophica (Bacterium, anoxygenic photosynthesis)
    
    indices = [2, 5, 6, 7, 9, 10, 12]
    
    # Convert each number to a string and join them with a comma
    print(",".join(map(str, indices)))

solve_task()