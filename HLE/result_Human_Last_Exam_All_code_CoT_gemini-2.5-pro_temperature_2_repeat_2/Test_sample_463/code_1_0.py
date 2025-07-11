def solve_photochemical_synthesis_query():
    """
    Identifies and prints the indices of taxa that perform photochemical synthesis.

    The decision for each taxon is based on established biological knowledge.
    "Photochemical synthesis" is interpreted as phototrophy: the use of light
    energy to drive metabolic processes like ATP synthesis or carbon fixation.
    This includes oxygenic photosynthesis, anoxygenic photosynthesis, and
    retinal-based phototrophy. It excludes animal processes like vision or
    Vitamin D synthesis, which are not for primary energy capture.
    """

    # A dictionary representing the taxa and their phototrophic capability.
    # True if they perform photochemical synthesis, False otherwise.
    taxa_synthesis_map = {
        1: False,  # Acanthella cavernosa (Sponge)
        2: True,   # Gloeochaete wittrockiana (Alga)
        3: False,  # Homo sapiens (Human)
        4: False,  # Riftia pachyptila (Tube worm, symbiosis ignored)
        5: True,   # Halapricum salinum (Archaea, bacteriorhodopsin)
        6: True,   # Aphanothece castagnei (Cyanobacterium)
        7: True,   # Baileya pleniradiata (Plant)
        8: False,  # Acanthella pulchra (Sponge)
        9: True,   # Ectothiorhodosinus mongolicus (Purple sulfur bacterium)
        10: True,  # Chlorobaculum tepidum (Green sulfur bacterium)
        11: False, # Stygichthys typhlops (Cavefish)
        12: True,  # Gemmatimonas phototrophica (Photosynthetic bacterium)
        13: False  # Myonera garretti (Deep-sea bivalve)
    }

    # Collect the indices of the taxa that are True
    phototrophic_indices = []
    for index, is_phototrophic in taxa_synthesis_map.items():
        if is_phototrophic:
            phototrophic_indices.append(str(index))
    
    # Join the indices with a comma for the final output
    result_string = ",".join(phototrophic_indices)
    
    print(result_string)

solve_photochemical_synthesis_query()