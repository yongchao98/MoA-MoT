def solve_task():
    """
    This function identifies the indices of taxa that undergo photochemical synthesis.
    The analysis is as follows:
    - 2 (Gloeochaete wittrockiana): Green alga, performs oxygenic photosynthesis.
    - 3 (Homo sapiens): Human, synthesizes Vitamin D using UV light.
    - 5 (Halapricum salinum): Archaea, uses bacteriorhodopsin for photophosphorylation (ATP synthesis).
    - 6 (Aphanothece castagnei): Cyanobacterium, performs oxygenic photosynthesis.
    - 7 (Baileya pleniradiata): Plant, performs oxygenic photosynthesis.
    - 9 (Ectothiorhodosinus mongolicus): Purple sulfur bacterium, performs anoxygenic photosynthesis.
    - 10 (Chlorobaculum tepidum): Green sulfur bacterium, performs anoxygenic photosynthesis.
    - 12 (Gemmatimonas phototrophica): Bacterium, performs anoxygenic photosynthesis.
    - The remaining taxa (1, 4, 8, 11, 13) are non-photosynthetic animals.
    """
    
    # List of indices for taxa that perform photochemical synthesis.
    photochemical_synthesizers_indices = [2, 3, 5, 6, 7, 9, 10, 12]
    
    # Convert numbers to strings to join them with a comma.
    result_string = ",".join(map(str, photochemical_synthesizers_indices))
    
    print(result_string)

solve_task()