def find_milkweed_mutualists():
    """
    Identifies the mutualists of Asclepias fascicularis from a predefined list.

    A mutualistic relationship is one where both organisms benefit. For Asclepias fascicularis
    (Narrowleaf Milkweed), this includes pollinators and defenders.

    The function uses a predefined dictionary that classifies the interaction of each organism
    with the milkweed.
    """

    # Data based on biological interactions:
    # Adults of bees, wasps, and butterflies are pollinators (mutualists).
    # Ants can be defenders, feeding on extrafloral nectaries (mutualists).
    # Larvae of butterflies are herbivores (not mutualists).
    # Larvae of bees and wasps are typically fed by adults and do not interact directly with the plant.
    interactions = {
        1: {'organism': 'Danaus plexipus (Monarch Butterfly)', 'stage': 'Adult', 'relationship': 'Mutualist'},
        2: {'organism': 'Megachile frigidus (Leafcutter Bee)', 'stage': 'Adult', 'relationship': 'Mutualist'},
        3: {'organism': 'Formica rufa (Red Wood Ant)', 'stage': 'Adult', 'relationship': 'Mutualist'},
        4: {'organism': 'Sphex ichneumoneus (Digger Wasp)', 'stage': 'Adult', 'relationship': 'Mutualist'},
        5: {'organism': 'Pepsis thisbe (Tarantula Hawk Wasp)', 'stage': 'Adult', 'relationship': 'Mutualist'},
        6: {'organism': 'Megachile ericetorum (Leafcutter Bee)', 'stage': 'Adult', 'relationship': 'Mutualist'},
        7: {'organism': 'Danaus plexipus (Monarch Butterfly)', 'stage': 'Larva', 'relationship': 'Herbivore'},
        8: {'organism': 'Megachile frigidus (Leafcutter Bee)', 'stage': 'Larva', 'relationship': 'No direct interaction'},
        9: {'organism': 'Formica rufa (Red Wood Ant)', 'stage': 'Larva', 'relationship': 'No direct interaction'},
        10: {'organism': 'Sphex ichneumoneus (Digger Wasp)', 'stage': 'Larva', 'relationship': 'No direct interaction'},
        11: {'organism': 'Pepsis thisbe (Tarantula Hawk Wasp)', 'stage': 'Larva', 'relationship': 'No direct interaction'},
        12: {'organism': 'Megachile ericetorum (Leafcutter Bee)', 'stage': 'Larva', 'relationship': 'No direct interaction'},
    }

    mutualist_indices = []
    for index, data in interactions.items():
        if data['relationship'] == 'Mutualist':
            mutualist_indices.append(str(index))

    if not mutualist_indices:
        print("none")
    else:
        result = ",".join(mutualist_indices)
        print(result)

find_milkweed_mutualists()