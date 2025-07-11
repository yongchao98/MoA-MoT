def find_photosynthetic_taxa():
    """
    This function identifies taxa that perform photochemical synthesis from a predefined list.

    The logic is based on biological knowledge:
    - 2: Gloeochaete wittrockiana is an alga (oxygenic photosynthesis).
    - 3: Homo sapiens synthesizes Vitamin D using sunlight.
    - 5: Halapricum salinum is a haloarchaeon, which can use light energy to synthesize ATP (photophosphorylation).
    - 6: Aphanothece castagnei is a cyanobacterium (oxygenic photosynthesis).
    - 7: Baileya pleniradiata is a plant (oxygenic photosynthesis).
    - 9: Ectothiorhodosinus mongolicus is a purple sulfur bacterium (anoxygenic photosynthesis).
    - 10: Chlorobaculum tepidum is a green sulfur bacterium (anoxygenic photosynthesis).
    - 12: Gemmatimonas phototrophica is a bacterium that performs anoxygenic photosynthesis.

    Other taxa are animals that do not perform photochemical synthesis themselves (ignoring symbionts).
    """
    
    # Indices of the taxa that undergo photochemical synthesis
    photochemical_indices = [2, 3, 5, 6, 7, 9, 10, 12]

    # Convert the list of numbers to a comma-separated string for printing
    result_string = ",".join(map(str, photochemical_indices))
    
    print(result_string)

find_photosynthetic_taxa()