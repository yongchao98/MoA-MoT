def find_photochemical_synthesizers():
    """
    Identifies and prints the indices of taxa that perform photochemical synthesis.

    The function stores a predefined list of taxa and their metabolic capabilities.
    It then filters this list to find organisms that use light energy to synthesize
    compounds as a normal part of their metabolism, ignoring symbiotic relationships.
    This includes various forms of photosynthesis and other specific photochemical reactions.
    """

    species_data = [
        # (Index, Name, Performs Photochemical Synthesis, Reason)
        (1, "Acanthella cavernosa", False, "Animal (sponge)"),
        (2, "Gloeochaete wittrockiana", True, "Alga (oxygenic photosynthesis)"),
        (3, "Homo sapiens", True, "Animal (Vitamin D synthesis via UV light)"),
        (4, "Riftia pachyptila", False, "Animal (chemosynthesis)"),
        (5, "Halapricum salinum", False, "Archaea (chemoorganotroph)"),
        (6, "Aphanothece castagnei", True, "Cyanobacterium (oxygenic photosynthesis)"),
        (7, "Baileya pleniradiata", True, "Plant (oxygenic photosynthesis)"),
        (8, "Acanthella pulchra", False, "Animal (sponge)"),
        (9, "Ectothiorhodosinus mongolicus", True, "Bacterium (anoxygenic photosynthesis)"),
        (10, "Chlorobaculum tepidum", True, "Bacterium (anoxygenic photosynthesis)"),
        (11, "Stygichthys typhlops", False, "Animal (cavefish)"),
        (12, "Gemmatimonas phototrophica", True, "Bacterium (anoxygenic photosynthesis)"),
        (13, "Myonera garretti", False, "Animal (bivalve)"),
    ]

    phototrophic_indices = []
    for index, name, is_phototrophic, reason in species_data:
        if is_phototrophic:
            phototrophic_indices.append(index)
    
    # Convert indices to strings and join with a comma
    result = ",".join(map(str, phototrophic_indices))
    print(result)

find_photochemical_synthesizers()