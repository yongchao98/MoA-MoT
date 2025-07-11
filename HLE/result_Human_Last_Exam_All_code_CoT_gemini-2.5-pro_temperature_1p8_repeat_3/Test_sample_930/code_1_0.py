def find_mutualists():
    """
    This function identifies the mutualists of Asclepias fascicularis from a predefined list.
    Mutualists are organisms that have a relationship where both species benefit.
    In this case, adult pollinators (butterflies, bees, wasps) and ant defenders are mutualists.
    Larvae that eat the plant (herbivores) or do not interact with it are not mutualists.
    """
    
    # Indices of organisms identified as mutualists
    # 1: Danaus plexippus (Adult) - Pollinator
    # 2: Megachile frigidus (Adult) - Pollinator
    # 3: Formica rufa (Adult) - Nectar feeder and defender
    # 4: Sphex ichneumoneus (Adult) - Pollinator
    # 5: Pepsis thisbe (Adult) - Pollinator
    # 6: Megachile ericetorum (Adult) - Pollinator
    mutualist_indices = [1, 2, 3, 4, 5, 6]
    
    # Convert indices to strings and join them with a comma
    result = ",".join(map(str, mutualist_indices))
    
    print(result)

find_mutualists()