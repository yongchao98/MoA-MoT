def find_mutualists():
    """
    Identifies the mutualists of Asclepias fascicularis from a predefined list.

    A mutualist is an organism that engages in a mutually beneficial relationship.
    In this context, mutualists are the adult insects that pollinate the plant
    in exchange for nectar.

    The predefined list is:
    Adults:
    1) Danaus plexipus (pollinator)
    2) Megachile frigidus (pollinator)
    3) Formica rufa (nectar robber, not a reliable pollinator)
    4) Sphex ichneumoneus (pollinator)
    5) Pepsis thisbe (pollinator)
    6) Megachile ericetorum (pollinator)

    Larvae:
    7) Danaus plexipus (herbivore, eats the plant)
    8) Megachile frigidus (develops in nest, no direct interaction)
    9) Formica rufa (develops in colony, no direct interaction)
    10) Sphex ichneumoneus (carnivore, no direct interaction)
    11) Pepsis thisbe (parasitoid, no direct interaction)
    12) Megachile ericetorum (develops in nest, no direct interaction)
    """

    # Indices of the organisms identified as mutualists
    mutualist_indices = [1, 2, 4, 5, 6]

    # Convert the list of numbers to a list of strings to join them
    mutualist_indices_str = [str(i) for i in mutualist_indices]

    # Join the indices with a comma
    result = ",".join(mutualist_indices_str)

    print(result)

find_mutualists()