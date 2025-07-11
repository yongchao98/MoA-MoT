def find_mutualists():
    """
    This function identifies the indices of organisms that have a mutualistic
    relationship with Asclepias fascicularis from a predefined list.

    Mutualism is a relationship where both species benefit. In this case:
    - Adult insects (butterflies, bees, wasps) act as pollinators while feeding on nectar.
    - Ants can act as protectors against herbivores while feeding on nectar.
    - Larvae are either herbivores (harming the plant) or do not interact with the plant at all.

    Therefore, the mutualists are all the adult forms listed.
    """
    # Indices of the mutualists based on biological analysis
    mutualist_indices = [1, 2, 3, 4, 5, 6]

    # Convert each index to a string for printing
    output_string = ",".join(map(str, mutualist_indices))

    print(output_string)

find_mutualists()