def find_species_with_reduced_pterostigmata():
    """
    Identifies dragonfly species with reduced pterostigmata based on their ecology.

    The key ecological trait associated with reduced pterostigmata in dragonflies is
    long-distance migration combined with a gliding/soaring flight style. Species
    that have mastered this flight (e.g., to cross oceans) have evolved wings
    optimized for soaring, which often includes a very small pterostigma.
    """

    # The full list of taxa provided by the user.
    all_taxa = [
        "Didymops transversa",
        "Urothemis edwarsi",
        "Macrodiplax balteata",
        "Pantala flavescens",
        "Orthetrum cancellatum",
        "Libelulla quadrimaculata",
        "Libelulla pulchela",
        "Sympetrum corruptum",
        "Celithemis elisa",
        "Tholymis tillarga"
    ]

    # Species known for being long-distance migratory "gliders" with reduced pterostigmata.
    # This list is based on established ecological and morphological knowledge.
    target_species = [
        "Macrodiplax balteata", # Known migrant and glider, often associated with Pantala.
        "Pantala flavescens",   # The quintessential example of a wandering glider.
        "Tholymis tillarga"     # A well-known migratory glider, often active at dusk.
    ]

    # Find the 1-based indices of the target species in the full list.
    indices = []
    for i, taxon in enumerate(all_taxa):
        if taxon in target_species:
            # The user's list is 1-based, so we add 1 to the 0-based index.
            indices.append(i + 1)
    
    # Sort the indices for a consistent output format.
    indices.sort()

    # Format the result as a comma-separated string.
    result_string = ",".join(map(str, indices))

    print(result_string)

find_species_with_reduced_pterostigmata()
<<<3,4,10>>>