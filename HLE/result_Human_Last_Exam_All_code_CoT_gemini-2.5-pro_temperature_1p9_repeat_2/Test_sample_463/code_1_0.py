def find_photochemical_synthesis_taxa():
    """
    Identifies which of a given list of taxa undergo photochemical synthesis
    as a primary metabolic process (phototrophy).
    """

    taxa = [
        {"index": 1, "name": "Acanthella cavernosa", "type": "Animal (Sponge)"},
        {"index": 2, "name": "Gloeochaete wittrockiana", "type": "Alga (Glaucophyte)"},
        {"index": 3, "name": "Homo sapiens", "type": "Animal (Mammal)"},
        {"index": 4, "name": "Riftia pachyptila", "type": "Animal (Tube Worm)"},
        {"index": 5, "name": "Halapricum salinum", "type": "Archaea (Halophile)"},
        {"index": 6, "name": "Aphanothece castagnei", "type": "Bacteria (Cyanobacterium)"},
        {"index": 7, "name": "Baileya pleniradiata", "type": "Plant"},
        {"index": 8, "name": "Acanthella pulchra", "type": "Animal (Sponge)"},
        {"index": 9, "name": "Ectothiorhodosinus mongolicus", "type": "Bacteria (Purple Sulfur)"},
        {"index": 10, "name": "Chlorobaculum tepidum", "type": "Bacteria (Green Sulfur)"},
        {"index": 11, "name": "Stygichthys typhlops", "type": "Animal (Cavefish)"},
        {"index": 12, "name": "Gemmatimonas phototrophica", "type": "Bacteria (Gemmatimonadetes)"},
        {"index": 13, "name": "Myonera garretti", "type": "Animal (Bivalve)"}
    ]

    # Keywords identifying groups known for phototrophy
    phototroph_keywords = [
        "Alga",
        "Plant",
        "Cyanobacterium",
        "Sulfur",         # Catches purple and green sulfur bacteria
        "Halophile",      # Catches archaea using retinal-based phototrophy
        "phototrophica"   # Catches the species named for being phototrophic
    ]

    result_indices = []

    # Iterate through each taxon and check if it's a phototroph
    for taxon in taxa:
        is_phototroph = False
        # Create a single string with the taxon's name and type for searching
        search_string = taxon["type"] + " " + taxon["name"]
        
        for keyword in phototroph_keywords:
            if keyword in search_string:
                is_phototroph = True
                break
        
        if is_phototroph:
            result_indices.append(str(taxon["index"]))

    if not result_indices:
        print("none")
    else:
        # Join the numbers into a comma-separated string for the final output
        output_string = ",".join(result_indices)
        print(output_string)

find_photochemical_synthesis_taxa()