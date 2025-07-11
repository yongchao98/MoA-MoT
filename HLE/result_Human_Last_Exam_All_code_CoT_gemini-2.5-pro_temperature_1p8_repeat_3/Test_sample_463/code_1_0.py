def solve_task():
    """
    Identifies which of the provided taxa undergo photochemical synthesis.
    """
    species = [
        "Acanthella cavernosa",
        "Gloeochaete wittrockiana",
        "Homo sapiens",
        "Riftia pachyptila",
        "Halapricum salinum",
        "Aphanothece castagnei",
        "Baileya pleniradiata",
        "Acanthella pulchra",
        "Ectothiorhodosinus mongolicus",
        "Chlorobaculum tepidum",
        "Stygichthys typhlops",
        "Gemmatimonas phototrophica",
        "Myonera garretti"
    ]

    # A boolean list indicating if the species at the corresponding index
    # performs photochemical synthesis.
    # 1. Sponge (Animal) -> No
    # 2. Algae -> Yes
    # 3. Human (Animal) -> No
    # 4. Tube Worm (Animal, Chemoautotroph) -> No
    # 5. Archaea (Phototroph) -> Yes
    # 6. Cyanobacteria -> Yes
    # 7. Plant -> Yes
    # 8. Sponge (Animal) -> No
    # 9. Purple Sulfur Bacteria -> Yes
    # 10. Green Sulfur Bacteria -> Yes
    # 11. Cavefish (Animal) -> No
    # 12. Gemmatimonadetes Bacteria -> Yes
    # 13. Bivalve (Animal) -> No
    performs_photochemical_synthesis = [
        False, True, False, False, True, True, True, False, True, True, False, True, False
    ]

    phototrophic_indices = []
    for i, is_phototrophic in enumerate(performs_photochemical_synthesis):
        if is_phototrophic:
            # The list of species is 1-indexed for the user
            phototrophic_indices.append(str(i + 1))
    
    if phototrophic_indices:
        result = ",".join(phototrophic_indices)
    else:
        result = "none"
        
    print(result)

solve_task()