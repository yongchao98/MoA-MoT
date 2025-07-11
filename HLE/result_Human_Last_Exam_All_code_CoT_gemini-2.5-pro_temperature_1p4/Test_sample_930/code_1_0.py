def find_milkweed_mutualists():
    """
    Identifies mutualists of Asclepias fascicularis from a predefined list.

    In the context of plant-insect interactions, mutualism typically involves
    pollination. Adult insects that feed on nectar and effectively pollinate the
    plant are considered mutualists. Herbivores and insects that do not interact
    with the plant in a beneficial way are excluded.
    """

    # Data representing each organism and its primary interaction with Asclepias fascicularis
    # 'pollinator': Mutualist (gets nectar, pollinates plant)
    # 'herbivore': Antagonist (eats plant)
    # 'nectar_robber': Not a mutualist (takes nectar, does not pollinate effectively)
    # 'no_direct_interaction': Larval stages that do not interact with the plant
    organisms = {
        1: {'species': 'Danaus plexipus', 'stage': 'Adult', 'interaction': 'pollinator'},
        2: {'species': 'Megachile frigidus', 'stage': 'Adult', 'interaction': 'pollinator'},
        3: {'species': 'Formica rufa', 'stage': 'Adult', 'interaction': 'nectar_robber'},
        4: {'species': 'Sphex ichneumoneus', 'stage': 'Adult', 'interaction': 'pollinator'},
        5: {'species': 'Pepsis thisbe', 'stage': 'Adult', 'interaction': 'pollinator'},
        6: {'species': 'Megachile ericetorum', 'stage': 'Adult', 'interaction': 'pollinator'},
        7: {'species': 'Danaus plexipus', 'stage': 'Larva', 'interaction': 'herbivore'},
        8: {'species': 'Megachile frigidus', 'stage': 'Larva', 'interaction': 'no_direct_interaction'},
        9: {'species': 'Formica rufa', 'stage': 'Larva', 'interaction': 'no_direct_interaction'},
        10: {'species': 'Sphex ichneumoneus', 'stage': 'Larva', 'interaction': 'no_direct_interaction'},
        11: {'species': 'Pepsis thisbe', 'stage': 'Larva', 'interaction': 'no_direct_interaction'},
        12: {'species': 'Megachile ericetorum', 'stage': 'Larva', 'interaction': 'no_direct_interaction'},
    }

    mutualist_indices = []
    for index, data in organisms.items():
        if data['interaction'] == 'pollinator':
            mutualist_indices.append(str(index))

    if mutualist_indices:
        print(",".join(mutualist_indices))
    else:
        print("none")

find_milkweed_mutualists()