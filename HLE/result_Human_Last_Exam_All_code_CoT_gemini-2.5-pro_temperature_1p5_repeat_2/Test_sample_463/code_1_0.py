def find_photosynthetic_taxa():
    """
    Identifies which taxa from a predefined list undergo photochemical synthesis.
    
    The function checks a list of organisms and their known metabolic characteristics.
    It identifies those that perform any form of photochemical synthesis,
    such as photosynthesis or bacteriorhodopsin-based phototrophy,
    as a normal part of their metabolism, ignoring symbiotic relationships.
    """
    
    taxa_data = {
        1: ("Acanthella cavernosa", False),  # Animal (sponge), heterotroph
        2: ("Gloeochaete wittrockiana", True),  # Algae, oxygenic photosynthesis
        3: ("Homo sapiens", True),  # Animal, Vitamin D synthesis via UV light
        4: ("Riftia pachyptila", False), # Animal (tube worm), chemoautotroph
        5: ("Halapricum salinum", True),  # Archaea, bacteriorhodopsin-based phototrophy
        6: ("Aphanothece castagnei", True), # Cyanobacteria, oxygenic photosynthesis
        7: ("Baileya pleniradiata", True),  # Plant, oxygenic photosynthesis
        8: ("Acanthella pulchra", False), # Animal (sponge), heterotroph
        9: ("Ectothiorhodosinus mongolicus", True), # Bacteria, anoxygenic photosynthesis
        10: ("Chlorobaculum tepidum", True), # Bacteria, anoxygenic photosynthesis
        11: ("Stygichthys typhlops", False), # Animal (cavefish), heterotroph
        12: ("Gemmatimonas phototrophica", True), # Bacteria, anoxygenic photosynthesis
        13: ("Myonera garretti", False)  # Animal (bivalve), heterotroph
    }
    
    photochemical_indices = []
    for index, (name, is_photochemical) in taxa_data.items():
        if is_photochemical:
            photochemical_indices.append(str(index))
            
    if not photochemical_indices:
        print("none")
    else:
        print(",".join(photochemical_indices))

find_photosynthetic_taxa()