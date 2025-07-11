def find_milkweed_mutualists():
    """
    Identifies mutualists of Asclepias fascicularis from a predefined list.

    A mutualistic relationship is one where both species benefit. For Asclepias,
    this typically means an insect acts as a pollinator (helping the plant reproduce)
    while it feeds on the plant's nectar (getting food).
    """

    # List of organisms with their index, name, and relationship to the plant.
    organisms = {
        1: ("Danaus plexipus (Adult)", "Pollinator", True),
        2: ("Megachile frigidus (Adult)", "Pollinator", True),
        3: ("Formica rufa (Adult)", "Nectar Robber", False),
        4: ("Sphex ichneumoneus (Adult)", "Pollinator", True),
        5: ("Pepsis thisbe (Adult)", "Pollinator", True),
        6: ("Megachile ericetorum (Adult)", "Pollinator", True),
        7: ("Danaus plexipus (Larvae)", "Herbivore (eats leaves)", False),
        8: ("Megachile frigidus (Larvae)", "No direct interaction", False),
        9: ("Formica rufa (Larvae)", "No direct interaction", False),
        10: ("Sphex ichneumoneus (Larvae)", "No direct interaction", False),
        11: ("Pepsis thisbe (Larvae)", "No direct interaction", False),
        12: ("Megachile ericetorum (Larvae)", "No direct interaction", False)
    }

    print("Step-by-step evaluation to find mutualists of Asclepias fascicularis:\n")
    print("A mutualist is an organism that provides a benefit to the plant (e.g., pollination) while receiving a benefit itself (e.g., nectar).")
    print("-" * 20)

    mutualist_indices = []
    
    # Iterate through the organisms and explain the reasoning.
    for index, (name, relationship, is_mutualist) in organisms.items():
        if is_mutualist:
            print(f"{index}) {name}: This adult insect is a pollinator. It gets nectar and helps the plant reproduce. -> Mutualist")
            mutualist_indices.append(index)
        else:
            print(f"{index}) {name}: This organism is an herbivore, a nectar robber, or has no direct interaction. -> Not a mutualist")

    # Sort indices for consistent output
    mutualist_indices.sort()

    print("\n" + "-" * 20)
    print("The final list of mutualists is composed of the adult pollinators.")
    # The "equation" is the collection of indices that satisfy the condition.
    print("The indices of all mutualists are:")
    
    # Final formatted output
    result_string = ",".join(map(str, mutualist_indices))
    print(result_string)

find_milkweed_mutualists()
<<<1,2,4,5,6>>>