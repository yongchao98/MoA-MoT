def find_mutualists():
    """
    This function identifies the indices of organisms that are mutualists
    with Asclepias fascicularis from a predefined list.

    A mutualistic relationship is one where both species benefit.
    In this context, mutualists are primarily pollinators or defenders.

    - Adult butterflies, bees, and wasps (1, 2, 4, 5, 6) are pollinators,
      getting nectar while helping the plant reproduce.
    - Ants (3) can act as defenders, protecting the plant from herbivores
      in exchange for nectar.
    - Larvae are either herbivores (7) or do not interact with the plant (8-12).
    """

    # List of indices corresponding to the mutualists
    mutualist_indices = [1, 2, 3, 4, 5, 6]

    # Convert the list of numbers to a list of strings for joining
    mutualist_indices_str = [str(i) for i in mutualist_indices]

    # Join the strings with a comma
    result = ",".join(mutualist_indices_str)

    print(result)

find_mutualists()