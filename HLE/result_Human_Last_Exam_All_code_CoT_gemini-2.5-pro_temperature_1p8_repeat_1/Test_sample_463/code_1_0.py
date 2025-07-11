def find_photosynthetic_taxa():
    """
    Identifies which of a given list of taxa undergo photochemical synthesis.
    """
    taxa_info = {
        1: ("Acanthella cavernosa", False, "Marine sponge, heterotroph"),
        2: ("Gloeochaete wittrockiana", True, "Green alga, performs oxygenic photosynthesis"),
        3: ("Homo sapiens", True, "Human, performs photochemical synthesis of Vitamin D"),
        4: ("Riftia pachyptila", False, "Tube worm, relies on chemosynthesis"),
        5: ("Halapricum salinum", True, "Halophilic archaeon, uses bacteriorhodopsin for phototrophy (ATP synthesis)"),
        6: ("Aphanothece castagnei", True, "Cyanobacterium, performs oxygenic photosynthesis"),
        7: ("Baileya pleniradiata", True, "Plant, performs oxygenic photosynthesis"),
        8: ("Acanthella pulchra", False, "Marine sponge, heterotroph"),
        9: ("Ectothiorhodosinus mongolicus", True, "Purple sulfur bacterium, performs anoxygenic photosynthesis"),
        10: ("Chlorobaculum tepidum", True, "Green sulfur bacterium, performs anoxygenic photosynthesis"),
        11: ("Stygichthys typhlops", False, "Cave fish, heterotroph"),
        12: ("Gemmatimonas phototrophica", True, "Bacterium, performs aerobic anoxygenic phototrophy"),
        13: ("Myonera garretti", False, "Bivalve mollusc, heterotroph"),
    }

    photochemical_indices = []
    for index, (name, performs_photochem, reason) in taxa_info.items():
        if performs_photochem:
            photochemical_indices.append(str(index))
            
    result = ",".join(photochemical_indices)
    print(result)

find_photosynthetic_taxa()