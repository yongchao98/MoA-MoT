def find_photochemical_taxa():
    """
    Identifies taxa that perform photochemical synthesis from a predefined list.

    The function stores a list of taxa with a boolean flag indicating whether they
    perform any photochemical process (e.g., photosynthesis, vitamin D synthesis,
    bacteriorhodopsin-based phototrophy) as part of their normal metabolism,
    ignoring symbiotic relationships. It then filters this list and prints the
    indices of the qualifying taxa in a comma-separated format.
    """
    taxa_data = {
        1: ("Acanthella cavernosa", False),  # Sponge (Animal), heterotroph
        2: ("Gloeochaete wittrockiana", True),   # Alga, performs photosynthesis
        3: ("Homo sapiens", True),            # Animal, performs Vitamin D synthesis via UV light
        4: ("Riftia pachyptila", False),      # Deep-sea worm (Animal), chemosynthesis is symbiotic
        5: ("Halapricum salinum", True),         # Archaea, performs photophosphorylation (bacteriorhodopsin)
        6: ("Aphanothece castagnei", True),   # Cyanobacterium, performs oxygenic photosynthesis
        7: ("Baileya pleniradiata", True),     # Plant, performs photosynthesis
        8: ("Acanthella pulchra", False),     # Sponge (Animal), heterotroph
        9: ("Ectothiorhodosinus mongolicus", True), # Bacterium, performs anoxygenic photosynthesis
        10: ("Chlorobaculum tepidum", True),  # Bacterium, performs anoxygenic photosynthesis
        11: ("Stygichthys typhlops", False),    # Cave fish (Animal), lives in darkness
        12: ("Gemmatimonas phototrophica", True), # Bacterium, performs anoxygenic photosynthesis
        13: ("Myonera garretti", False)         # Bivalve (Animal), heterotroph
    }

    photochemical_indices = []
    for index, (name, performs_synthesis) in taxa_data.items():
        if performs_synthesis:
            photochemical_indices.append(str(index))

    # The problem asks for the species indices separated by ",".
    result_string = ",".join(photochemical_indices)
    print(result_string)

find_photochemical_taxa()