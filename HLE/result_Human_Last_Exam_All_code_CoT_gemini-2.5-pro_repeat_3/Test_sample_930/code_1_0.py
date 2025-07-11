def find_mutualists():
    """
    Identifies mutualists of Asclepias fascicularis from a predefined list.

    A mutualistic relationship is one where both species benefit. In this case,
    it involves insects that pollinate the plant or protect it in exchange for food (nectar).
    - Adult butterflies, bees, and wasps are pollinators.
    - Adult ants can be protectors.
    - Larvae that eat the plant are herbivores (harmful).
    - Larvae that don't interact with the plant are neutral.
    """

    # Data representing each organism and its relationship with Asclepias fascicularis
    # 'relationship' can be 'mutualism', 'herbivory' (harmful), or 'none' (neutral)
    organisms = {
        1: {'name': 'Danaus plexipus (Adult)', 'relationship': 'mutualism'}, # Pollinator
        2: {'name': 'Megachile frigidus (Adult)', 'relationship': 'mutualism'}, # Pollinator
        3: {'name': 'Formica rufa (Adult)', 'relationship': 'mutualism'}, # Nectar-feeding protector
        4: {'name': 'Sphex ichneumoneus (Adult)', 'relationship': 'mutualism'}, # Pollinator
        5: {'name': 'Pepsis thisbe (Adult)', 'relationship': 'mutualism'}, # Pollinator
        6: {'name': 'Megachile ericetorum (Adult)', 'relationship': 'mutualism'}, # Pollinator
        7: {'name': 'Danaus plexipus (Larva)', 'relationship': 'herbivory'}, # Feeds on leaves
        8: {'name': 'Megachile frigidus (Larva)', 'relationship': 'none'}, # In nest, no interaction
        9: {'name': 'Formica rufa (Larva)', 'relationship': 'none'}, # In nest, no interaction
        10: {'name': 'Sphex ichneumoneus (Larva)', 'relationship': 'none'}, # In burrow, no interaction
        11: {'name': 'Pepsis thisbe (Larva)', 'relationship': 'none'}, # In burrow, no interaction
        12: {'name': 'Megachile ericetorum (Larva)', 'relationship': 'none'} # In nest, no interaction
    }

    mutualist_indices = []
    for index, data in organisms.items():
        if data['relationship'] == 'mutualism':
            mutualist_indices.append(str(index))

    if mutualist_indices:
        # The final list of indices for the mutualists
        result = ",".join(mutualist_indices)
        print(result)
    else:
        print("none")

find_mutualists()