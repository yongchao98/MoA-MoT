def find_mutualists():
    """
    This function identifies and prints the indices of organisms that are
    mutualists of Asclepias fascicularis from the provided list.

    A mutualistic relationship is one where both species benefit. In this case,
    the mutualists are the adult insects that pollinate the milkweed flowers
    in exchange for nectar.
    """
    # List of indices of mutualists based on biological knowledge.
    # Adults (1-6): Butterflies, bees, and wasps are pollinators (mutualists). Ants are nectar robbers.
    # Larvae (7-12): The monarch larva is a herbivore (antagonistic). The other larvae do not interact with the plant.
    mutualist_indices = [1, 2, 4, 5, 6]

    # The format required is a string of indices separated by commas.
    # We convert each integer index to a string before joining.
    output_string = ",".join(map(str, mutualist_indices))

    print(output_string)

find_mutualists()