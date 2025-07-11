def find_species_with_reduced_pterostigmata():
    """
    Identifies dragonfly species with reduced pterostigmata based on their ecology.

    The primary ecological indicator for a reduced pterostigma in dragonflies
    is a highly aerial, long-distance migratory lifestyle that involves extensive
    gliding. This wing modification reduces profile drag, enhancing flight efficiency.
    Species that are primarily perchers or non-migratory "cruisers" typically have
    more prominent pterostigmata for flight stability.
    """
    
    # The full list of taxa provided by the user.
    all_taxa = [
        "Didymops transversa",      # Stream Cruiser: A cruiser, not a long-distance migrant.
        "Urothemis edwarsi",        # Elegant Basker: A percher, not a noted migrant.
        "Macrodiplax balteata",     # Marl Pennant: A known long-distance migrant.
        "Pantala flavescens",       # Globe Skimmer: The world's most widespread and famous migratory dragonfly.
        "Orthetrum cancellatum",    # Black-tailed Skimmer: A percher, not a primary migrant.
        "Libelulla quadrimaculata", # Four-spotted Chaser: Can be migratory but has a prominent pterostigma.
        "Libelulla pulchela",       # Twelve-spotted Skimmer: A territorial percher.
        "Sympetrum corruptum",      # Variegated Meadowhawk: A well-known migratory species.
        "Celithemis elisa",         # Calico Pennant: A percher, not a migrant.
        "Tholymis tillarga"         # Coral-tailed Cloudwing: A strong crepuscular migrant.
    ]

    # Species with a migratory/gliding ecology indicating reduced pterostigmata.
    migratory_gliders = {
        "Macrodiplax balteata",
        "Pantala flavescens",
        "Sympetrum corruptum",
        "Tholymis tillarga"
    }

    found_indices = []
    # Enumerate through the list to get both the 1-based index and the name.
    for index, taxon in enumerate(all_taxa, 1):
        if taxon in migratory_gliders:
            found_indices.append(index)

    # Format the output as requested.
    if not found_indices:
        print("none")
    else:
        # Convert the list of integer indices to a comma-separated string.
        result_string = ",".join(map(str, sorted(found_indices)))
        print(result_string)

# Execute the function to find and print the indices.
find_species_with_reduced_pterostigmata()