def find_species_with_reduced_pterostigmata():
    """
    Identifies dragonfly species expected to have reduced pterostigmata based on their ecology.

    The pterostigma adds weight to the wing tip and is involved in controlling wing torsion
    and flutter during flapping flight. In species that are extreme long-distance migrants
    and specialize in a "gliding" or "sailing" flight style, minimizing weight at the
    wingtips is a key adaptation for energy efficiency. Therefore, these species are
    expected to have a reduced pterostigma.

    This function analyzes a list of species and identifies those known for this
    migratory, gliding ecology.
    """
    # A list of dictionaries representing the species and their ecological flight type.
    # The 'is_migratory_glider' flag is set for species in the Pantalinae subfamily,
    # which are famous for their gliding flight and extensive migrations.
    species_data = [
        {"name": "Didymops transversa", "is_migratory_glider": False},       # Percher, not a long-distance migrant.
        {"name": "Urothemis edwarsi", "is_migratory_glider": False},         # Percher.
        {"name": "Macrodiplax balteata", "is_migratory_glider": False},      # Migratory, but not a classic "glider".
        {"name": "Pantala flavescens", "is_migratory_glider": True},        # The quintessential "Wandering Glider", extreme migrant.
        {"name": "Orthetrum cancellatum", "is_migratory_glider": False},     # Percher.
        {"name": "Libelulla quadrimaculata", "is_migratory_glider": False}, # Can be migratory, but primarily a percher.
        {"name": "Libelulla pulchela", "is_migratory_glider": False},       # Percher.
        {"name": "Sympetrum corruptum", "is_migratory_glider": False},       # Migratory, but not specialized as a "glider".
        {"name": "Celithemis elisa", "is_migratory_glider": False},          # Percher.
        {"name": "Tholymis tillarga", "is_migratory_glider": True}          # "Coral-tailed Cloudwing", a crepuscular glider and migrant.
    ]

    reduced_pterostigmata_indices = []
    for i, species in enumerate(species_data):
        # The list is 1-based in the prompt, so we use i + 1 for the index.
        if species["is_migratory_glider"]:
            reduced_pterostigmata_indices.append(str(i + 1))

    # Print the indices separated by a comma, as requested.
    print(",".join(reduced_pterostigmata_indices))

find_species_with_reduced_pterostigmata()