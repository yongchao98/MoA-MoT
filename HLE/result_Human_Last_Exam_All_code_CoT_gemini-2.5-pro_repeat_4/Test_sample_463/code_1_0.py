def solve_task():
    """
    Identifies which of a given list of taxa perform photochemical synthesis.

    The function maintains a list of the taxa with a boolean flag indicating
    whether they are known to be phototrophic (i.e., use light energy for synthesis).
    It then filters this list to find the indices of the phototrophic organisms
    and prints them as a comma-separated string.
    """
    taxa_data = [
        (1, "Acanthella cavernosa", False),             # Animal (sponge)
        (2, "Gloeochaete wittrockiana", True),          # Alga (photosynthetic)
        (3, "Homo sapiens", False),                      # Animal (heterotroph)
        (4, "Riftia pachyptila", False),                # Animal (relies on chemosynthesis)
        (5, "Halapricum salinum", True),               # Archaea (photoheterotroph)
        (6, "Aphanothece castagnei", True),            # Cyanobacteria (photosynthetic)
        (7, "Baileya pleniradiata", True),             # Plant (photosynthetic)
        (8, "Acanthella pulchra", False),               # Animal (sponge)
        (9, "Ectothiorhodosinus mongolicus", True),   # Bacteria (anoxygenic photosynthesis)
        (10, "Chlorobaculum tepidum", True),           # Bacteria (anoxygenic photosynthesis)
        (11, "Stygichthys typhlops", False),           # Animal (fish)
        (12, "Gemmatimonas phototrophica", True),      # Bacteria (anoxygenic photosynthesis)
        (13, "Myonera garretti", False)                # Animal (bivalve)
    ]

    # Collect the indices of all phototrophic organisms
    phototrophic_indices = []
    for index, name, is_phototrophic in taxa_data:
        if is_phototrophic:
            phototrophic_indices.append(str(index))

    # Join the indices with a comma and print the result
    result_string = ",".join(phototrophic_indices)
    print(result_string)

solve_task()