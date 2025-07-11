def find_mutualists():
    """
    Identifies the mutualists of Asclepias fascicularis from a predefined list.

    In ecology, mutualism is a relationship where both species benefit. For the
    plant Asclepias fascicularis (milkweed), mutualists are primarily pollinators
    that feed on nectar and help the plant reproduce.

    - Danaus plexippus (adults, 1), Megachile frigidus (adults, 2),
      Sphex ichneumoneus (adults, 4), Pepsis thisbe (adults, 5), and
      Megachile ericetorum (adults, 6) are all known pollinators of milkweed.
      They gain nectar (food) and the plant gets pollinated. This is mutualism.
    - Danaus plexippus larvae (7) are herbivores that eat the plant's leaves, which is harmful
      to the plant, not mutualistic.
    - The other larvae (8, 9, 10, 11, 12) do not interact with the plant.
    - Formica rufa (ants, 3) are typically nectar robbers, taking nectar without providing a
      significant pollination service, so they are not considered mutualists.
    """
    mutualist_indices = [1, 2, 4, 5, 6]

    # The print statement below will format the output as requested.
    # It converts each number in the list to a string and then joins them
    # with a comma and a space.
    print(",".join(map(str, mutualist_indices)))

find_mutualists()