import sys

def find_species_with_reduced_pterostigmata():
    """
    Identifies dragonfly species expected to have reduced pterostigmata based on their ecology.

    The primary ecological factor correlated with reduced pterostigmata is a highly
    migratory, "glider" lifestyle. This flight style, common in species that travel
    long distances, often over open water, favors aerodynamic efficiency where the
    stabilizing weight of a large pterostigma may be less critical or adapted differently.

    This script codifies this biological principle by filtering a list of species
    for those known to fit this ecological profile.
    """

    # A data store containing the names and relevant ecological information for each species.
    # The 'ecology' field focuses on flight style and migratory behavior.
    species_database = [
        {'id': 1, 'name': 'Didymops transversa', 'ecology': 'Riverine cruiser, not a specialized glider.'},
        {'id': 2, 'name': 'Urothemis edwarsi', 'ecology': 'Resident of still waters, not typically a long-distance migrant or glider.'},
        {'id': 3, 'name': 'Macrodiplax balteata', 'ecology': 'Highly migratory and often found over open water; considered a glider.'},
        {'id': 4, 'name': 'Pantala flavescens', 'ecology': 'The quintessential "Globe Skimmer," a trans-oceanic migrant and specialized glider.'},
        {'id': 5, 'name': 'Orthetrum cancellatum', 'ecology': 'Territorial skimmer, not a specialized long-distance migrant or glider.'},
        {'id': 6, 'name': 'Libelulla quadrimaculata', 'ecology': 'Known for migratory flights, but not classified as a specialized glider.'},
        {'id': 7, 'name': 'Libelulla pulchela', 'ecology': 'A common perching skimmer, not a specialized glider.'},
        {'id': 8, 'name': 'Sympetrum corruptum', 'ecology': 'A migratory species, but not primarily known as a glider.'},
        {'id': 9, 'name': 'Celithemis elisa', 'ecology': 'A "percher" species, not a glider.'},
        {'id': 10, 'name': 'Tholymis tillarga', 'ecology': 'A highly migratory, crepuscular species belonging to the "glider" (Trameinae) group.'}
    ]

    candidate_indices = []
    # Iterate through the database to find species that match the "glider" profile.
    for species in species_database:
        # We search for the keyword 'glider' as the key identifier for the relevant ecology.
        if 'glider' in species['ecology'].lower():
            candidate_indices.append(str(species['id']))

    # Format the output as a comma-separated string of indices, or 'none'.
    if candidate_indices:
        print(",".join(candidate_indices))
    else:
        print("none")

find_species_with_reduced_pterostigmata()