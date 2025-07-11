def solve():
    """
    Identifies which of a given list of taxa undergo photochemical synthesis.

    The function holds a predefined list of taxa and their known metabolic capabilities
    regarding the use of light energy. It filters this list to find the organisms
    that perform photochemical synthesis (e.g., photosynthesis, retinal-based phototrophy)
    and prints their corresponding indices in a comma-separated format.
    """
    
    # List of tuples: (index, name, performs_photochemical_synthesis)
    taxa_info = [
        (1, "Acanthella cavernosa", False),           # Animal (Sponge)
        (2, "Gloeochaete wittrockiana", True),    # Algae (Oxygenic photosynthesis)
        (3, "Homo sapiens", False),                # Animal (Human)
        (4, "Riftia pachyptila", False),           # Animal (Tube worm, relies on chemosynthesis)
        (5, "Halapricum salinum", True),          # Archaea (Retinal-based phototrophy)
        (6, "Aphanothece castagnei", True),       # Bacteria (Cyanobacteria, oxygenic photosynthesis)
        (7, "Baileya pleniradiata", True),        # Plant (Oxygenic photosynthesis)
        (8, "Acanthella pulchra", False),          # Animal (Sponge)
        (9, "Ectothiorhodosinus mongolicus", True), # Bacteria (Purple sulfur, anoxygenic photosynthesis)
        (10, "Chlorobaculum tepidum", True),      # Bacteria (Green sulfur, anoxygenic photosynthesis)
        (11, "Stygichthys typhlops", False),       # Animal (Cavefish)
        (12, "Gemmatimonas phototrophica", True),  # Bacteria (Anoxygenic photosynthesis)
        (13, "Myonera garretti", False)            # Animal (Bivalve)
    ]

    # Find the indices of taxa that perform photochemical synthesis
    phototrophic_indices = []
    for index, name, is_phototrophic in taxa_info:
        if is_phototrophic:
            phototrophic_indices.append(str(index))

    # Join the indices with a comma and print the result
    result_string = ",".join(phototrophic_indices)
    print(result_string)

solve()