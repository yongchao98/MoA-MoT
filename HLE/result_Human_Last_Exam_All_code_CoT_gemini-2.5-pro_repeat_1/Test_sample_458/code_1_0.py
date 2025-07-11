def find_species_with_reduced_pterostigmata():
    """
    Identifies dragonfly species from a predefined list that are expected
    to have reduced pterostigmata based on their ecology.

    The primary ecological factor associated with reduced pterostigmata in
    dragonflies is a migratory, "glider" lifestyle. This flight style is
    common in species that undertake long-distance, often trans-oceanic,
    migrations. The reduction in the pterostigma is an adaptation for
    more efficient soaring and gliding flight.

    Based on this, the following species are selected:
    - 3) Macrodiplax balteata: A known migrant and strong flier, often dispersing over coastal waters.
    - 4) Pantala flavescens: The "Globe Skimmer," famous for being the most widespread and migratory dragonfly in the world.
    - 10) Tholymis tillarga: A "glider" known for its migratory and crepuscular (dawn/dusk) flight habits.
    """
    
    # Indices of the species with the glider ecology
    species_indices = [3, 4, 10]
    
    # Sort the indices for consistent output
    species_indices.sort()
    
    # Format the result as a comma-separated string
    result = ",".join(map(str, species_indices))
    
    print(result)

find_species_with_reduced_pterostigmata()
<<<3,4,10>>>