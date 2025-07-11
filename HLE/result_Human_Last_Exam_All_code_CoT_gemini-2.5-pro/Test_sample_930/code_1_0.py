def find_mutualists():
    """
    This function identifies the indices of mutualists of Asclepias fascicularis
    from a predefined list and prints them in the required format.
    """
    # Adults:
    # 1) Danaus plexipus (Monarch Butterfly) - Pollinator -> Mutualist
    # 2) Megachile frigidus (Leafcutter Bee) - Pollinator -> Mutualist
    # 3) Formica rufa (Red Wood Ant) - Protector/Nectar feeder -> Mutualist
    # 4) Sphex ichneumoneus (Great Golden Digger Wasp) - Pollinator -> Mutualist
    # 5) Pepsis thisbe (Tarantula Hawk Wasp) - Pollinator -> Mutualist
    # 6) Megachile ericetorum (Leafcutter Bee) - Pollinator -> Mutualist
    
    # Larvae:
    # 7) Danaus plexipus (Monarch Caterpillar) - Herbivore -> Not a mutualist
    # 8-12) Other larvae do not interact with the plant in a mutualistic way.

    mutualist_indices = [1, 2, 3, 4, 5, 6]
    
    # Convert numbers to strings and join with a comma
    result = ",".join(map(str, mutualist_indices))
    
    print(result)

find_mutualists()