def find_mutualists():
    """
    Identifies the mutualists of Asclepias fascicularis from a predefined list.

    In ecology, mutualism is an interaction between two species where both benefit.
    For flowering plants like Asclepias fascicularis (milkweed), this typically
    involves pollination: the insect gets food (nectar), and the plant gets pollinated.

    - Adult butterflies (1), bees (2, 6), and nectar-feeding wasps (4, 5) are pollinators and thus mutualists.
    - Ants (3) are generally considered nectar robbers, not effective pollinators.
    - The Monarch larva (7) is an herbivore; it eats the plant.
    - The other larvae (8, 9, 10, 11, 12) do not directly interact with the plant in a mutualistic way.
    """
    
    # List of indices corresponding to the mutualists
    mutualist_indices = [1, 2, 4, 5, 6]
    
    # The final answer should be the indices separated by a comma.
    # We convert each number to a string and then join them with a comma.
    output = ",".join(map(str, mutualist_indices))
    
    print(output)

find_mutualists()