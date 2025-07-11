def find_milkweed_mutualists():
    """
    Identifies mutualists of Asclepias fascicularis from a predefined list.

    This function encodes biological knowledge about the interactions between
    various insects and milkweed. A mutualist is defined here as an organism
    that provides a benefit to the plant, such as pollination or defense,
    while receiving a benefit, like nectar.
    """
    
    # Step 1 & 2: Define the organisms and their interactions.
    # Interaction types: 'pollinator', 'defender', 'herbivore', 'parasitoid', 'no_direct_interaction'
    interactions = {
        1:  ("Danaus plexipus", "Adult", "pollinator"),
        2:  ("Megachile frigidus", "Adult", "pollinator"),
        3:  ("Formica rufa", "Adult", "defender"),
        4:  ("Sphex ichneumoneus", "Adult", "pollinator"),
        5:  ("Pepsis thisbe", "Adult", "pollinator"),
        6:  ("Megachile ericetorum", "Adult", "pollinator"),
        7:  ("Danaus plexipus", "Larvae", "herbivore"),
        8:  ("Megachile frigidus", "Larvae", "no_direct_interaction"),
        9:  ("Formica rufa", "Larvae", "no_direct_interaction"),
        10: ("Sphex ichneumoneus", "Larvae", "parasitoid"),
        11: ("Pepsis thisbe", "Larvae", "parasitoid"),
        12: ("Megachile ericetorum", "Larvae", "no_direct_interaction"),
    }

    # Step 3: Filter for mutualistic interactions.
    mutualistic_interactions = {'pollinator', 'defender'}
    mutualist_indices = []

    for index, data in interactions.items():
        interaction_type = data[2]
        if interaction_type in mutualistic_interactions:
            mutualist_indices.append(str(index))

    # Step 4: Format and print the final result.
    if mutualist_indices:
        result = ",".join(mutualist_indices)
        print(result)
    else:
        print("none")

find_milkweed_mutualists()