def find_mutualists():
    """
    This function identifies and prints the indices of mutualists of Asclepias fascicularis
    from a predefined list.
    
    A mutualistic relationship is identified for:
    - Adult insects that act as pollinators (receiving nectar).
    - Adult insects that act as protectors (receiving nectar).
    
    A non-mutualistic relationship is identified for:
    - Larvae that are herbivores (e.g., Monarch caterpillars).
    - Larvae that do not interact with the plant at all.
    """
    
    # List of indices for species that are mutualists
    # 1: Danaus plexipus (adult) - Pollinator
    # 2: Megachile frigidus (adult) - Pollinator
    # 3: Formica rufa (adult) - Protector
    # 4: Sphex ichneumoneus (adult) - Pollinator
    # 5: Pepsis thisbe (adult) - Pollinator
    # 6: Megachile ericetorum (adult) - Pollinator
    mutualist_indices = [1, 2, 3, 4, 5, 6]
    
    # Convert the list of numbers to a list of strings for joining
    mutualist_indices_str = [str(i) for i in mutualist_indices]
    
    # Join the strings with a comma and print the result
    result = ",".join(mutualist_indices_str)
    print(result)

find_mutualists()