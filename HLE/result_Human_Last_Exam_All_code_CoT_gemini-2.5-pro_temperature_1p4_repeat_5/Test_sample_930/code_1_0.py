def find_mutualists():
    """
    Identifies mutualists of Asclepias fascicularis from a predefined list.

    A mutualistic relationship for a plant involves receiving a benefit,
    typically pollination. The analysis considers the organism's life stage,
    its interaction with the plant, and its native geographic range.

    - Asclepias fascicularis is native to North America.
    - Adult pollinators that feed on nectar provide a benefit.
    - Larvae are typically herbivores (harmful) or develop elsewhere and do
      not interact with the plant.
    - Species not native to North America do not have a co-evolved relationship.
    """
    organisms = [
        # 1: Adult monarch butterfly - Pollinator, correct range. MUTUALIST.
        {'index': 1, 'name': 'Danaus plexipus', 'stage': 'Adult', 'is_mutualist': True},
        # 2: Adult leaf-cutter bee - Pollinator, correct range. MUTUALIST.
        {'index': 2, 'name': 'Megachile frigidus', 'stage': 'Adult', 'is_mutualist': True},
        # 3: Adult red wood ant - Not native to North America. NOT MUTUALIST.
        {'index': 3, 'name': 'Formica rufa', 'stage': 'Adult', 'is_mutualist': False},
        # 4: Adult great golden digger wasp - Pollinator, correct range. MUTUALIST.
        {'index': 4, 'name': 'Sphex ichneumoneus', 'stage': 'Adult', 'is_mutualist': True},
        # 5: Adult tarantula hawk wasp - Pollinator, correct range. MUTUALIST.
        {'index': 5, 'name': 'Pepsis thisbe', 'stage': 'Adult', 'is_mutualist': True},
        # 6: Adult leaf-cutter bee - Not native to North America. NOT MUTUALIST.
        {'index': 6, 'name': 'Megachile ericetorum', 'stage': 'Adult', 'is_mutualist': False},
        # 7: Monarch larva - Herbivore, eats the plant. NOT MUTUALIST.
        {'index': 7, 'name': 'Danaus plexipus', 'stage': 'Larva', 'is_mutualist': False},
        # 8: Bee larva - Develops in a nest, no interaction. NOT MUTUALIST.
        {'index': 8, 'name': 'Megachile frigidus', 'stage': 'Larva', 'is_mutualist': False},
        # 9: Ant larva - Develops in a nest, no interaction. NOT MUTUALIST.
        {'index': 9, 'name': 'Formica rufa', 'stage': 'Larva', 'is_mutualist': False},
        # 10: Wasp larva - Carnivore, develops in a burrow. NOT MUTUALIST.
        {'index': 10, 'name': 'Sphex ichneumoneus', 'stage': 'Larva', 'is_mutualist': False},
        # 11: Wasp larva - Parasitoid on spiders. NOT MUTUALIST.
        {'index': 11, 'name': 'Pepsis thisbe', 'stage': 'Larva', 'is_mutualist': False},
        # 12: Bee larva - Develops in a nest, no interaction. NOT MUTUALIST.
        {'index': 12, 'name': 'Megachile ericetorum', 'stage': 'Larva', 'is_mutualist': False},
    ]

    mutualist_indices = [str(org['index']) for org in organisms if org['is_mutualist']]

    if mutualist_indices:
        result = ",".join(mutualist_indices)
    else:
        result = "none"

    print(result)

find_mutualists()
<<<1,2,4,5>>>