def find_mutualists():
    """
    This function identifies the indices of mutualists of Asclepias fascicularis
    from a predefined list and prints them.

    The list of potential mutualists is:
    Adults:
    1) Danaus plexipus (Pollinator)
    2) Megachile frigidus (Pollinator)
    3) Formica rufa (Defender/Nectar-drinker)
    4) Sphex ichneumoneus (Pollinator)
    5) Pepsis thisbe (Pollinator)
    6) Megachile ericetorum (Pollinator)
    Larvae:
    7) Danaus plexipus (Herbivore)
    8) Megachile frigidus (In nest, no direct interaction)
    9) Formica rufa (In nest, no direct interaction)
    10) Sphex ichneumoneus (Parasitoid, in nest)
    11) Pepsis thisbe (Parasitoid, on host)
    12) Megachile ericetorum (In nest, no direct interaction)

    A mutualist relationship is one where both species benefit. For the plant,
    this usually involves pollination or defense in exchange for nectar.
    """

    # Indices of the organisms that are mutualists
    mutualist_indices = [1, 2, 3, 4, 5, 6]

    # Convert the list of numbers to a list of strings
    mutualist_strings = [str(index) for index in mutualist_indices]

    # Join the strings with a comma
    result = ",".join(mutualist_strings)

    # Print the final result
    print(result)

find_mutualists()