def find_species_with_reduced_pterostigmata():
    """
    Identifies dragonfly species from a predefined list that are expected
    to have reduced pterostigmata based on their ecology.

    The primary ecological trait associated with reduced pterostigmata is
    long-distance, wind-assisted migration (a "glider" lifestyle).
    """
    taxa_list = [
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

    # These species are known long-distance "gliders" and have reduced pterostigmata.
    species_with_reduced_pterostigmata = [
        "Macrodiplax balteata",
        "Pantala flavescens",
        "Tholymis tillarga"
    ]

    # Find the 1-based indices of the selected species in the original list
    indices = []
    for species in species_with_reduced_pterostigmata:
        if species in taxa_list:
            # The question uses 1-based indexing
            indices.append(taxa_list.index(species) + 1)
    
    indices.sort()
    
    # Format the indices as a comma-separated string
    result_string = ",".join(map(str, indices))
    
    print(result_string)

find_species_with_reduced_pterostigmata()