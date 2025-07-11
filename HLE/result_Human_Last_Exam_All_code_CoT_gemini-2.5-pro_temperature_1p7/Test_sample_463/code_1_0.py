def find_photosynthetic_taxa():
    """
    Identifies and prints the indices of taxa that perform photochemical synthesis.
    """
    # A list of tuples, where each tuple contains the species name and a boolean
    # indicating if it performs any form of photochemical synthesis.
    taxa_data = [
        ("Acanthella cavernosa", False),        # 1. Sponge (Animal)
        ("Gloeochaete wittrockiana", True),     # 2. Algae (Photosynthesis)
        ("Homo sapiens", True),                 # 3. Human (Vitamin D synthesis)
        ("Riftia pachyptila", False),           # 4. Tube worm (Chemosynthesis)
        ("Halapricum salinum", True),           # 5. Archaea (Bacteriorhodopsin phototrophy)
        ("Aphanothece castagnei", True),        # 6. Cyanobacteria (Photosynthesis)
        ("Baileya pleniradiata", True),         # 7. Plant (Photosynthesis)
        ("Acanthella pulchra", False),          # 8. Sponge (Animal)
        ("Ectothiorhodosinus mongolicus", True),# 9. Bacteria (Anoxygenic Photosynthesis)
        ("Chlorobaculum tepidum", True),        # 10. Bacteria (Anoxygenic Photosynthesis)
        ("Stygichthys typhlops", False),        # 11. Cavefish (Animal)
        ("Gemmatimonas phototrophica", True),   # 12. Bacteria (Anoxygenic Photosynthesis)
        ("Myonera garretti", False)             # 13. Bivalve (Animal)
    ]

    photosynthetic_indices = []
    # Iterate through the list with an index starting from 1
    for index, (taxon, is_photosynthetic) in enumerate(taxa_data, 1):
        if is_photosynthetic:
            photosynthetic_indices.append(str(index))

    # Join the indices into a comma-separated string and print
    result = ",".join(photosynthetic_indices)
    print(result)

find_photosynthetic_taxa()