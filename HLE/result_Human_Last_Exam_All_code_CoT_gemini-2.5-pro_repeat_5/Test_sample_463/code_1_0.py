def find_photochemical_synthesizers():
    """
    Identifies which of the given taxa undergo photochemical synthesis.
    """
    species_data = {
        1: ("Acanthella cavernosa", False),
        2: ("Gloeochaete wittrockiana", True),
        3: ("Homo sapiens", False),
        4: ("Riftia pachyptila", False),
        5: ("Halapricum salinum", False),
        6: ("Aphanothece castagnei", True),
        7: ("Baileya pleniradiata", True),
        8: ("Acanthella pulchra", False),
        9: ("Ectothiorhodosinus mongolicus", True),
        10: ("Chlorobaculum tepidum", True),
        11: ("Stygichthys typhlops", False),
        12: ("Gemmatimonas phototrophica", True),
        13: ("Myonera garretti", False)
    }

    phototrophic_indices = []
    for index, (name, is_phototrophic) in species_data.items():
        if is_phototrophic:
            phototrophic_indices.append(str(index))

    if not phototrophic_indices:
        print("none")
    else:
        result = ",".join(phototrophic_indices)
        print(result)

find_photochemical_synthesizers()