def find_species_with_reduced_pterostigmata():
    """
    Identifies dragonfly species from a predefined list that are expected to have
    reduced pterostigmata based on their ecology.

    The primary ecological trait linked to reduced pterostigmata is a lifestyle
    of long-distance migration and extensive gliding/soaring. The pterostigma,
    a weight at the wing tip, prevents flutter in powered flight but is less
    necessary and adds weight for gliders.

    Based on established entomological knowledge:
    - Macrodiplax balteata (3): A strong coastal migrant known for gliding.
    - Pantala flavescens (4): The 'Globe Skimmer', the world's most widespread
      and longest-distance migratory dragonfly, a quintessential glider.
    - Tholymis tillarga (10): A pan-tropical migrant also known for its gliding flight.
    """
    
    # The indices of the species from the user's list that fit the criteria.
    species_indices = [3, 4, 10]
    
    # Convert each integer index to a string for the final output.
    # The .join() method concatenates them with a comma separator.
    output_string = ",".join(map(str, species_indices))
    
    print(output_string)

find_species_with_reduced_pterostigmata()