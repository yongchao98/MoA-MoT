def find_glider_species():
    """
    Identifies dragonfly species from a predefined list that are expected to have
    reduced pterostigmata based on their ecology.

    The species with this trait are typically long-distance migratory "gliders",
    a lifestyle characteristic of the tribe Trameini. This group has evolved a
    suite of features for efficient gliding, including broad hindwings and
    relatively small pterostigmata.

    The identified species from the list are:
    - 3: Macrodiplax balteata
    - 4: Pantala flavescens
    - 10: Tholymis tillarga
    """
    
    # Indices of the species with a "glider" ecology.
    glider_indices = [3, 4, 10]
    
    # Convert the list of integers to a list of strings for joining
    # This ensures the output format is exactly as requested.
    output_string = ",".join(map(str, glider_indices))
    
    print(output_string)

find_glider_species()