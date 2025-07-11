def find_milkweed_mutualists():
    """
    Identifies mutualists of Asclepias fascicularis from a predefined list.

    In ecology, mutualism is a relationship where both species benefit. For a
    flowering plant like milkweed, the primary mutualism with insects is pollination.
    This code evaluates a list of insects to determine which ones are pollinators.
    """

    # Data representing the organisms and their relationship with milkweed.
    # 'pollinator' is a mutualistic relationship.
    # 'herbivore' and 'nectar_robber' are antagonistic.
    # 'indirect' means the stage does not directly interact with the plant.
    organisms = [
        {'index': 1, 'name': 'Danaus plexipus', 'stage': 'Adult', 'relationship': 'pollinator'},
        {'index': 2, 'name': 'Megachile frigidus', 'stage': 'Adult', 'relationship': 'pollinator'},
        {'index': 3, 'name': 'Formica rufa', 'stage': 'Adult', 'relationship': 'nectar_robber'},
        {'index': 4, 'name': 'Sphex ichneumoneus', 'stage': 'Adult', 'relationship': 'pollinator'},
        {'index': 5, 'name': 'Pepsis thisbe', 'stage': 'Adult', 'relationship': 'pollinator'},
        {'index': 6, 'name': 'Megachile ericetorum', 'stage': 'Adult', 'relationship': 'pollinator'},
        {'index': 7, 'name': 'Danaus plexipus', 'stage': 'Larva', 'relationship': 'herbivore'},
        {'index': 8, 'name': 'Megachile frigidus', 'stage': 'Larva', 'relationship': 'indirect'},
        {'index': 9, 'name': 'Formica rufa', 'stage': 'Larva', 'relationship': 'indirect'},
        {'index': 10, 'name': 'Sphex ichneumoneus', 'stage': 'Larva', 'relationship': 'indirect'},
        {'index': 11, 'name': 'Pepsis thisbe', 'stage': 'Larva', 'relationship': 'indirect'},
        {'index': 12, 'name': 'Megachile ericetorum', 'stage': 'Larva', 'relationship': 'indirect'}
    ]

    mutualist_indices = []
    for org in organisms:
        if org['relationship'] == 'pollinator':
            mutualist_indices.append(str(org['index']))

    if mutualist_indices:
        print(",".join(mutualist_indices))
    else:
        print("none")

find_milkweed_mutualists()