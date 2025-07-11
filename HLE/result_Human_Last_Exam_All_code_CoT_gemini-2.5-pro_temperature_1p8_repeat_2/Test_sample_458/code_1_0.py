def find_species_with_reduced_pterostigmata():
    """
    Identifies dragonfly species with reduced pterostigmata based on their ecology.

    The rationale is as follows:
    - Pterostigmata are essential for preventing wing flutter in strong, gliding fliers.
    - Long-distance migrants like Pantala flavescens (4), Macrodiplax balteata (3), and Sympetrum corruptum (8) require well-developed pterostigmata.
    - Strong territorial fliers like Didymops transversa (1) and Libellula species (6, 7) also have prominent pterostigmata.
    - Celithemis elisa (9): Belongs to the "pennants," known for a weaker, fluttering flight, which correlates with smaller pterostigmata.
    - Tholymis tillarga (10): This is a crepuscular (active at dawn/dusk) species. This ecology is linked to a flight style that does not require large pterostigmata, which are characteristically reduced in this genus.

    Therefore, the species expected to have reduced pterostigmata are Celithemis elisa and Tholymis tillarga.
    """
    
    # Indices of the identified species
    indices = [9, 10]
    
    # Format the output string by converting numbers to strings and joining with a comma
    output_string = ",".join(map(str, indices))
    
    print(output_string)

find_species_with_reduced_pterostigmata()