def find_milkweed_mutualists():
    """
    This function identifies and prints the indices of mutualists for Asclepias fascicularis
    based on their ecological roles.

    In this context, a mutualist is defined as a pollinator.
    - Adult butterflies, bees, and wasps that feed on nectar are pollinators.
    - Larvae that eat the plant are herbivores (antagonists).
    - Larvae that develop elsewhere have no interaction.
    - Ants are typically considered nectar robbers or minor pollinators at best, not primary mutualists.
    """
    # A dictionary representing the relationship of each organism with the milkweed.
    organism_roles = {
        1: "Pollinator",  # Danaus plexipus (Adult) - Feeds on nectar, pollinates
        2: "Pollinator",  # Megachile frigidus (Adult) - Feeds on nectar/pollen, pollinates
        3: "Antagonist",  # Formica rufa (Adult) - Nectar robber
        4: "Pollinator",  # Sphex ichneumoneus (Adult) - Feeds on nectar, pollinates
        5: "Pollinator",  # Pepsis thisbe (Adult) - Feeds on nectar, pollinates
        6: "Pollinator",  # Megachile ericetorum (Adult) - Feeds on nectar/pollen, pollinates
        7: "Herbivore",   # Danaus plexipus (Larva) - Eats leaves
        8: "No interaction", # Megachile frigidus (Larva) - Develops in a nest
        9: "No interaction", # Formica rufa (Larva) - Develops in a nest
        10: "No interaction",# Sphex ichneumoneus (Larva) - Develops on a host
        11: "No interaction",# Pepsis thisbe (Larva) - Develops on a host
        12: "No interaction" # Megachile ericetorum (Larva) - Develops in a nest
    }

    mutualist_indices = []
    for index, role in organism_roles.items():
        if role == "Pollinator":
            mutualist_indices.append(str(index))

    if mutualist_indices:
        result_string = ", ".join(mutualist_indices)
        # Printing the final list of numbers as requested.
        print(result_string)
    else:
        print("none")

find_milkweed_mutualists()