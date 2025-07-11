def find_photochemical_synthesizers():
    """
    Identifies and prints the indices of taxa that perform photochemical synthesis.

    The list of taxa and their metabolic capabilities is pre-determined based on biological knowledge.
    The function filters this list and prints the indices of the qualifying organisms.
    """
    # Each tuple: (index, scientific_name, performs_photochemical_synthesis)
    taxa_data = [
        (1, "Acanthella cavernosa", False),
        (2, "Gloeochaete wittrockiana", True),
        (3, "Homo sapiens", True),
        (4, "Riftia pachyptila", False),
        (5, "Halapricum salinum", True),
        (6, "Aphanothece castagnei", True),
        (7, "Baileya pleniradiata", True),
        (8, "Acanthella pulchra", False),
        (9, "Ectothiorhodosinus mongolicus", True),
        (10, "Chlorobaculum tepidum", True),
        (11, "Stygichthys typhlops", False),
        (12, "Gemmatimonas phototrophica", True),
        (13, "Myonera garretti", False)
    ]

    phototroph_indices = []
    for index, name, is_phototroph in taxa_data:
        if is_phototroph:
            phototroph_indices.append(str(index))

    # The final list of indices, formatted as a comma-separated string.
    # The numbers in the output are: 2, 3, 5, 6, 7, 9, 10, 12
    result = ",".join(phototroph_indices)
    
    print(result)

if __name__ == "__main__":
    find_photochemical_synthesizers()