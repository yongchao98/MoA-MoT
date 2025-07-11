def find_photochemical_taxa():
    """
    Identifies and prints the indices of taxa that undergo photochemical synthesis
    from a predefined list. The analysis is based on established biological knowledge.
    """
    # A dictionary mapping the index to the taxon name and a boolean indicating
    # if it performs any photochemical synthesis (e.g., photosynthesis, Vitamin D synthesis).
    # Symbiotic relationships are ignored as per the instructions.
    taxa_data = {
        1: ("Acanthella cavernosa", False),           # Animal (Sponge), heterotroph
        2: ("Gloeochaete wittrockiana", True),      # Algae, performs oxygenic photosynthesis
        3: ("Homo sapiens", True),                  # Animal, performs photochemical Vitamin D synthesis
        4: ("Riftia pachyptila", False),            # Animal (Tube worm), chemosynthesis-based, no light
        5: ("Halapricum salinum", False),           # Archaea, genome shows no phototrophy genes
        6: ("Aphanothece castagnei", True),         # Cyanobacteria, performs oxygenic photosynthesis
        7: ("Baileya pleniradiata", True),          # Plant, performs oxygenic photosynthesis
        8: ("Acanthella pulchra", False),           # Animal (Sponge), heterotroph
        9: ("Ectothiorhodosinus mongolicus", True), # Bacteria, performs anoxygenic photosynthesis
        10: ("Chlorobaculum tepidum", True),        # Bacteria, performs anoxygenic photosynthesis
        11: ("Stygichthys typhlops", False),        # Animal (Cavefish), heterotroph, no light
        12: ("Gemmatimonas phototrophica", True),   # Bacteria, performs anoxygenic photosynthesis
        13: ("Myonera garretti", False)             # Animal (Clam), heterotroph
    }

    photochemical_indices = []
    for index, data in taxa_data.items():
        is_photochemical = data[1]
        if is_photochemical:
            photochemical_indices.append(str(index))

    if photochemical_indices:
        print(",".join(photochemical_indices))
    else:
        print("none")

find_photochemical_taxa()