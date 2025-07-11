def find_milkweed_mutualists():
    """
    Identifies mutualists of Asclepias fascicularis from a given list.

    A mutualistic relationship is one where both species benefit. For Asclepias, this
    includes:
    - Pollinators (adult insects): They get nectar and the plant gets pollinated.
    - Defenders (e.g., ants): They get nectar and protect the plant from herbivores.

    Relationships that are NOT mutualistic include:
    - Herbivory (e.g., larvae eating leaves): The insect benefits, but the plant is harmed.
    - No direct interaction: The larval stages of many of these insects do not interact
      with the plant at all.
    """

    # Data representing each organism, its life stage, and its primary interaction with milkweed.
    organisms = {
        1: {'name': 'Danaus plexipus (Adult)', 'interaction': 'pollinator'},
        2: {'name': 'Megachile frigidus (Adult)', 'interaction': 'pollinator'},
        3: {'name': 'Formica rufa (Adult)', 'interaction': 'defender'},
        4: {'name': 'Sphex ichneumoneus (Adult)', 'interaction': 'pollinator'},
        5: {'name': 'Pepsis thisbe (Adult)', 'interaction': 'pollinator'},
        6: {'name': 'Megachile ericetorum (Adult)', 'interaction': 'pollinator'},
        7: {'name': 'Danaus plexipus (Larva)', 'interaction': 'herbivore'},
        8: {'name': 'Megachile frigidus (Larva)', 'interaction': 'none'},
        9: {'name': 'Formica rufa (Larva)', 'interaction': 'none'},
        10: {'name': 'Sphex ichneumoneus (Larva)', 'interaction': 'none'},
        11: {'name': 'Pepsis thisbe (Larva)', 'interaction': 'none'},
        12: {'name': 'Megachile ericetorum (Larva)', 'interaction': 'none'}
    }

    mutualist_indices = []
    for index, data in organisms.items():
        # Check if the interaction is beneficial to the plant (pollination or defense).
        if data['interaction'] in ['pollinator', 'defender']:
            mutualist_indices.append(str(index))

    if not mutualist_indices:
        print("none")
    else:
        # Format the output as a comma-separated string
        result = ",".join(mutualist_indices)
        print(result)

find_milkweed_mutualists()