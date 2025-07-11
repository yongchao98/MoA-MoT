def identify_migratory_gliders():
    """
    Identifies dragonfly species from a list that are expected to have reduced
    pterostigmata due to their "glider" or "wind-rider" ecology.

    The pterostigma is a wing cell that adds mass to the wingtip, helping to
    prevent flight-destabilizing vibrations. However, in highly specialized
    long-distance migratory dragonflies, known as "gliders," other aerodynamic
    adaptations, such as a broadened hindwing base for soaring, are more
    prominent. This specialized flight style can correlate with relatively
    reduced pterostigmata.

    The species from the list that fit this ecological profile are:
    - Macrodiplax balteata (3): A long-distance coastal migrant.
    - Pantala flavescens (4): The "Globe Skimmer," the most famous trans-oceanic migrant.
    - Tholymis tillarga (10): A crepuscular migrant that shares a similar gliding ecology.
    """
    
    species_list = [
        "Didymops transversa",       # 1
        "Urothemis edwarsi",         # 2
        "Macrodiplax balteata",      # 3
        "Pantala flavescens",        # 4
        "Orthetrum cancellatum",     # 5
        "Libelulla quadrimaculata",  # 6
        "Libelulla pulchela",        # 7
        "Sympetrum corruptum",       # 8
        "Celithemis elisa",          # 9
        "Tholymis tillarga"          # 10
    ]

    # List of species known for a "glider" lifestyle
    glider_species = [
        "Macrodiplax balteata",
        "Pantala flavescens",
        "Tholymis tillarga"
    ]

    # Find the 1-based indices of the identified species
    indices = []
    for i, species_name in enumerate(species_list):
        if species_name in glider_species:
            indices.append(str(i + 1))
    
    # Format and print the result as a comma-separated string
    result = ",".join(indices)
    print(result)

identify_migratory_gliders()