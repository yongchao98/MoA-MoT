def find_mutualists():
    """
    This function identifies the mutualists of Asclepias fascicularis from a predefined list.

    In the relationship with a plant like Asclepias fascicularis (milkweed), mutualists
    are organisms that engage in a mutually beneficial interaction. Typically, these are
    adult pollinators who consume nectar for food and, in the process, transfer pollen,
    aiding in the plant's reproduction.

    - Danaus plexippus (Monarch Butterfly) adults (1) are well-known pollinators.
    - Megachile species (leafcutter bees) adults (2, 6) are effective pollinators.
    - Sphex ichneumoneus (digger wasp) adults (4) are pollinators.
    - Pepsis thisbe (tarantula hawk wasp) adults (5) are also pollinators.
    - Formica rufa (ants) (3) are often nectar robbers and not effective pollinators.
    - The larvae of all listed insects are not mutualists. The Monarch larva (7) is a
      herbivore, and the other larvae (8-12) do not directly interact with the plant.

    Therefore, the indices corresponding to mutualists are identified.
    """
    mutualist_indices = [1, 2, 4, 5, 6]

    # The print statement formats the output as a comma-separated string
    # as requested by the user prompt.
    print(",".join(map(str, mutualist_indices)))

find_mutualists()