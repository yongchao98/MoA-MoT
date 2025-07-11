def find_migratory_dragonflies():
    """
    Identifies dragonflies likely to have reduced pterostigmata based on migratory ecology.
    
    The following species are known long-distance migrants, a lifestyle associated with
    reduced pterostigmata for improved gliding efficiency:
    
    3: Macrodiplax balteata (Wandering Pennant)
    4: Pantala flavescens (Wandering Glider)
    8: Sympetrum corruptum (Variegated Meadowhawk)
    10: Tholymis tillarga (Coral-tailed Cloudwing)
    """
    
    # Indices of the species with ecologies that suggest reduced pterostigmata.
    species_indices = [3, 4, 8, 10]
    
    # Format the output as a comma-separated string.
    result = ",".join(map(str, species_indices))
    
    print(result)

find_migratory_dragonflies()