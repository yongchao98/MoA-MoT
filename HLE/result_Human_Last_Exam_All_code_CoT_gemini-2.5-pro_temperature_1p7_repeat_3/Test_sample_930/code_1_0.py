def find_mutualists():
    """
    This function identifies and prints the indices of mutualists of Asclepias fascicularis from a predefined list.
    
    A mutualistic relationship is where both organisms benefit. For Asclepias, this is primarily pollination.
    
    - Mutualists are the adult insects that act as pollinators while feeding on nectar:
        1) Danaus plexipus (Monarch butterfly adult)
        2) Megachile frigidus (Leafcutter bee adult)
        4) Sphex ichneumoneus (Digger wasp adult)
        5) Pepsis thisbe (Tarantula hawk wasp adult)
        6) Megachile ericetorum (Leafcutter bee adult)
        
    - Not mutualists are:
        3) Formica rufa (Ants are poor pollinators and may protect plant pests).
        7) Danaus plexipus larvae (Herbivores that eat the plant).
        8-12) All other larvae (Do not interact with the plant directly).
    """
    
    # List of indices of the mutualists
    mutualist_indices = [1, 2, 4, 5, 6]
    
    # Convert the list of integers to a list of strings to use the join method
    mutualist_indices_str = [str(i) for i in mutualist_indices]
    
    # Join the elements with a comma to create the final string
    result = ",".join(mutualist_indices_str)
    
    print("The indices of all mutualists are:")
    print(result)

find_mutualists()