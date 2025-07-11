def find_species_with_reduced_pterostigmata():
    """
    Identifies dragonfly species from a predefined list that are expected
    to have reduced pterostigmata based on their ecology.

    The primary factor for reduced pterostigmata is a long-distance, migratory,
    "glider" lifestyle. Species that are exceptionally strong flyers and soar over
    vast distances often exhibit this adaptation.

    Based on established entomological knowledge:
    - Macrodiplax balteata (Index 3): A well-known migratory species, often found far out at sea.
    - Pantala flavescens (Index 4): The "Wandering Glider" or "Globe Skimmer," this is the most widespread and migratory dragonfly, famous for its gliding flight and reduced pterostigmata.
    - Tholymis tillarga (Index 10): Another pantropical "glider" known for its migratory habits, often flying with Pantala.
    
    The other species are generally considered "perchers," "chasers," or non-migratory "cruisers"
    and would be expected to have more prominent pterostigmata.
    """
    
    species_indices = [3, 4, 10]
    
    # Format the result as a comma-separated string
    result_string = ",".join(map(str, species_indices))
    
    print(f"The indices of the species expected to have reduced pterostigmata are: {result_string}")

find_species_with_reduced_pterostigmata()