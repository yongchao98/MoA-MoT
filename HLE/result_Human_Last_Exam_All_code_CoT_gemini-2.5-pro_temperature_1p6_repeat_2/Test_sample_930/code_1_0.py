def find_milkweed_mutualists():
    """
    This function identifies the mutualists of Asclepias fascicularis from a predefined list.

    A mutualist is an organism that has a mutually beneficial relationship with another.
    In this case, adult insects that pollinate the milkweed while feeding on its nectar are mutualists.
    The plant gets pollinated, and the insect gets food.
    
    The following indices correspond to mutualists:
    - 1: Danaus plexipus (Adult) - Pollinator
    - 2: Megachile frigidus (Adult) - Pollinator
    - 4: Sphex ichneumoneus (Adult) - Pollinator
    - 5: Pepsis thisbe (Adult) - Pollinator
    - 6: Megachile ericetorum (Adult) - Pollinator

    Larvae that eat the plant (e.g., Danaus plexipus larva) are herbivores, not mutualists.
    Other listed larvae and ants are not considered mutualists in this context.
    """
    
    mutualist_indices = [1, 2, 4, 5, 6]

    # The final answer is the list of indices, separated by a comma.
    # We will construct the output string to show each number explicitly.
    output_string = ",".join(map(str, mutualist_indices))
    
    print(output_string)

find_milkweed_mutualists()