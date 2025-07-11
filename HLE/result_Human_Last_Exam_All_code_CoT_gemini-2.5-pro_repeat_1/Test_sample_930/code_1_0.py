def find_mutualists():
    """
    Identifies the indices of organisms that are mutualists of Asclepias fascicularis
    from a predefined list.

    A mutualistic relationship in this context involves the insect acting as a pollinator
    in exchange for a food resource like nectar.

    - Adult butterflies (1), bees (2, 6), and nectar-feeding wasps (4, 5) are pollinators.
    - Ants (3) are nectar robbers and not effective pollinators.
    - Monarch larvae (7) are herbivores, which harms the plant.
    - All other larvae (8-12) do not directly interact with the plant in a beneficial way.
    """
    
    # Indices of the mutualists based on biological analysis
    mutualist_indices = [1, 2, 4, 5, 6]
    
    # Convert the list of numbers to a comma-separated string for output
    # The map(str, ...) function converts each number in the list to a string
    # The ",".join(...) function concatenates these strings with a comma in between
    output_string = ",".join(map(str, mutualist_indices))
    
    print(output_string)

find_mutualists()