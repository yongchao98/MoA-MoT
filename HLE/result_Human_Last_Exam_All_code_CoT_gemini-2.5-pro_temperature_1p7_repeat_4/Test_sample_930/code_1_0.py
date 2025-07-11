def find_mutualists():
    """
    Identifies and prints the indices of mutualists of Asclepias fascicularis
    from a predefined list.
    
    The list of organisms is:
    Adults:
    1) Danaus plexipus (pollinator)
    2) Megachile frigidus (pollinator)
    3) Formica rufa (defender/nectar feeder)
    4) Sphex ichneumoneus (pollinator)
    5) Pepsis thisbe (pollinator)
    6) Megachile ericetorum (pollinator)
    Larvae:
    7) Danaus plexipus (herbivore)
    8) Megachile frigidus (no direct interaction)
    9) Formica rufa (no direct interaction)
    10) Sphex ichneumoneus (no direct interaction)
    11) Pepsis thisbe (no direct interaction)
    12) Megachile ericetorum (no direct interaction)
    """
    # All listed adults are mutualists (pollinators or defenders).
    # None of the listed larvae are mutualists.
    mutualist_indices = [1, 2, 3, 4, 5, 6]
    
    # Convert the list of integers to a list of strings for joining
    indices_as_strings = [str(i) for i in mutualist_indices]
    
    # Join the strings with a comma and print the result
    print(",".join(indices_as_strings))

find_mutualists()