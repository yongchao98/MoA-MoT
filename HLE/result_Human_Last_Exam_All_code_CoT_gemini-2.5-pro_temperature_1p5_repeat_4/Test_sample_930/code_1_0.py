def solve_mutualism():
    """
    This function identifies the mutualists of Asclepias fascicularis from a predefined list.

    A mutualist is an organism that has a mutually beneficial relationship with another.
    In this case, we analyze each organism's interaction with the milkweed plant.

    - 1) Danaus plexipus (Adult): Pollinator. Mutualist.
    - 2) Megachile frigidus (Adult): Pollinator. Mutualist.
    - 3) Formica rufa (Adult): Nectar feeder and plant protector. Mutualist.
    - 4) Sphex ichneumoneus (Adult): Pollinator. Mutualist.
    - 5) Pepsis thisbe (Adult): Pollinator. Mutualist.
    - 6) Megachile ericetorum (Adult): Pollinator. Mutualist.
    - 7) Danaus plexipus (Larva): Herbivore (eats the plant). Not a mutualist.
    - 8) Megachile frigidus (Larva): Develops in a nest, fed by the adult. Not a mutualist.
    - 9) Formica rufa (Larva): Develops in a nest. Not a mutualist.
    - 10) Sphex ichneumoneus (Larva): Parasitoid on other insects. Not a mutualist.
    - 11) Pepsis thisbe (Larva): Parasitoid on other insects. Not a mutualist.
    - 12) Megachile ericetorum (Larva): Develops in a nest. Not a mutualist.

    The final list of mutualist indices is collected and printed.
    """
    mutualist_indices = [1, 2, 3, 4, 5, 6]

    # Check if the list is empty
    if not mutualist_indices:
        print("none")
    else:
        # Convert each integer index to a string for joining
        result_string = ",".join(map(str, mutualist_indices))
        print(result_string)

solve_mutualism()