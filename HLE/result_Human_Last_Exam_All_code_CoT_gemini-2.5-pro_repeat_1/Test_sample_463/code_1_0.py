def solve_photochemical_synthesis_query():
    """
    Identifies which of the given taxa perform photochemical synthesis.
    The list contains the index, name, and a boolean indicating if they perform
    any form of photochemical synthesis (e.g., photosynthesis, retinal-based phototrophy, Vitamin D synthesis).
    """
    taxa = [
        (1, "Acanthella cavernosa", False),          # Sponge (heterotroph)
        (2, "Gloeochaete wittrockiana", True),       # Algae (photosynthesis)
        (3, "Homo sapiens", True),                  # Human (Vitamin D synthesis)
        (4, "Riftia pachyptila", False),             # Tube worm (chemosynthesis via symbionts)
        (5, "Halapricum salinum", True),             # Archaea (bacteriorhodopsin phototrophy)
        (6, "Aphanothece castagnei", True),          # Cyanobacteria (photosynthesis)
        (7, "Baileya pleniradiata", True),           # Plant (photosynthesis)
        (8, "Acanthella pulchra", False),            # Sponge (heterotroph)
        (9, "Ectothiorhodosinus mongolicus", True),  # Purple sulfur bacteria (anoxygenic photosynthesis)
        (10, "Chlorobaculum tepidum", True),         # Green sulfur bacteria (anoxygenic photosynthesis)
        (11, "Stygichthys typhlops", False),         # Cave fish (heterotroph)
        (12, "Gemmatimonas phototrophica", True),    # Bacterium (anoxygenic photosynthesis)
        (13, "Myonera garretti", False)              # Bivalve (heterotroph)
    ]

    phototrophic_indices = []
    for index, name, is_phototrophic in taxa:
        if is_phototrophic:
            phototrophic_indices.append(str(index))
    
    # The final 'equation' is the comma-separated list of indices.
    result = ",".join(phototrophic_indices)
    print(result)

solve_photochemical_synthesis_query()