def solve_photochemical_synthesis_task():
    """
    Identifies which of the given taxa undergo any process of photochemical synthesis
    as part of their normal metabolic processes, ignoring symbiotic relationships.
    """

    taxa_data = {
        1: ("Acanthella cavernosa", False, "Sponge, heterotrophic animal."),
        2: ("Gloeochaete wittrockiana", True, "Green alga, performs oxygenic photosynthesis."),
        3: ("Homo sapiens", True, "Human, performs photochemical synthesis of Vitamin D in the skin."),
        4: ("Riftia pachyptila", False, "Tube worm, lives in deep sea without light, relies on chemosynthesis."),
        5: ("Halapricum salinum", True, "Archaea, performs retinal-based phototrophy (using light to pump protons)."),
        6: ("Aphanothece castagnei", True, "Cyanobacterium, performs oxygenic photosynthesis."),
        7: ("Baileya pleniradiata", True, "Plant, performs oxygenic photosynthesis."),
        8: ("Acanthella pulchra", False, "Sponge, heterotrophic animal."),
        9: ("Ectothiorhodosinus mongolicus", True, "Purple sulfur bacterium, performs anoxygenic photosynthesis."),
        10: ("Chlorobaculum tepidum", True, "Green sulfur bacterium, performs anoxygenic photosynthesis."),
        11: ("Stygichthys typhlops", False, "Cavefish, lives in darkness, heterotrophic animal."),
        12: ("Gemmatimonas phototrophica", True, "Bacterium, performs anoxygenic photosynthesis."),
        13: ("Myonera garretti", False, "Bivalve, heterotrophic animal.")
    }

    positive_indices = []
    for index, (name, performs_synthesis, reason) in taxa_data.items():
        if performs_synthesis:
            positive_indices.append(str(index))

    # The final output should be the indices separated by ","
    result = ",".join(positive_indices)
    print(result)

solve_photochemical_synthesis_task()
<<<2,3,5,6,7,9,10,12>>>