def find_mutualists():
    """
    Identifies the mutualists of Asclepias fascicularis from a predefined list.

    A mutualistic relationship is identified for:
    - Adult insects that act as pollinators (receiving nectar).
    - Adult insects that act as protectors (receiving nectar).

    Interactions that are not mutualistic:
    - Larval herbivory (antagonistic).
    - Life stages that do not interact with the plant.
    """
    # Indices of the mutualists based on biological interactions
    # 1: Danaus plexipus (adult) - Pollinator
    # 2: Megachile frigidus (adult) - Pollinator
    # 3: Formica rufa (adult) - Protector/Nectivore
    # 4: Sphex ichneumoneus (adult) - Pollinator
    # 5: Pepsis thisbe (adult) - Pollinator
    # 6: Megachile ericetorum (adult) - Pollinator
    mutualist_indices = [1, 2, 3, 4, 5, 6]

    # Convert the list of integers to a list of strings
    mutualist_indices_str = [str(i) for i in mutualist_indices]

    # Join the string elements with a comma
    result = ",".join(mutualist_indices_str)

    print(result)

find_mutualists()