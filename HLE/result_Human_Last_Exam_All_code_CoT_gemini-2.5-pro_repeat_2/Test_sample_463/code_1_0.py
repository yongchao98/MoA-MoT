def find_photochemical_synthesizers():
    """
    Identifies taxa that perform photochemical synthesis from a predefined list.

    The function iterates through a list of taxa, each with a flag indicating
    if it performs photochemical synthesis based on biological knowledge.
    It collects the indices of the qualifying taxa and prints them as a
    comma-separated string.
    """
    taxa_data = [
        {"index": 1, "name": "Acanthella cavernosa", "photochemical_synthesis": False},
        {"index": 2, "name": "Gloeochaete wittrockiana", "photochemical_synthesis": True},
        {"index": 3, "name": "Homo sapiens", "photochemical_synthesis": True},
        {"index": 4, "name": "Riftia pachyptila", "photochemical_synthesis": False},
        {"index": 5, "name": "Halapricum salinum", "photochemical_synthesis": True},
        {"index": 6, "name": "Aphanothece castagnei", "photochemical_synthesis": True},
        {"index": 7, "name": "Baileya pleniradiata", "photochemical_synthesis": True},
        {"index": 8, "name": "Acanthella pulchra", "photochemical_synthesis": False},
        {"index": 9, "name": "Ectothiorhodosinus mongolicus", "photochemical_synthesis": True},
        {"index": 10, "name": "Chlorobaculum tepidum", "photochemical_synthesis": True},
        {"index": 11, "name": "Stygichthys typhlops", "photochemical_synthesis": False},
        {"index": 12, "name": "Gemmatimonas phototrophica", "photochemical_synthesis": True},
        {"index": 13, "name": "Myonera garretti", "photochemical_synthesis": False}
    ]

    phototrophic_indices = []
    for taxon in taxa_data:
        if taxon["photochemical_synthesis"]:
            phototrophic_indices.append(str(taxon["index"]))

    if not phototrophic_indices:
        result = "none"
    else:
        result = ",".join(phototrophic_indices)
    
    print(result)

if __name__ == "__main__":
    find_photochemical_synthesizers()