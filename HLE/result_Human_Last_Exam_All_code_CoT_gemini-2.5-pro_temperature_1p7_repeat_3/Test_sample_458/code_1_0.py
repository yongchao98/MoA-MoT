def find_species_with_reduced_pterostigmata():
    """
    Identifies dragonfly species from a predefined list expected to have 
    reduced pterostigmata based on their ecology.

    The key ecological trait associated with reduced pterostigmata is a 
    highly specialized, long-distance, trans-oceanic migratory lifestyle. 
    This flight behavior favors wings adapted for maximum gliding efficiency.

    The analysis identifies the following species:
    - 3) Macrodiplax balteata: A known trans-oceanic migrant.
    - 4) Pantala flavescens: The "Wandering Glider," famous for the longest insect migrations.
    - 10) Tholymis tillarga: A widespread long-distance migrant, often flying over the sea.
    """
    
    # The indices of the species identified through ecological analysis.
    species_indices = [3, 4, 10]

    # Convert the list of integer indices to a list of strings
    species_indices_str = [str(index) for index in species_indices]

    # Join the string indices with a comma
    final_answer = ",".join(species_indices_str)

    # Print the final answer as requested
    print(final_answer)

find_species_with_reduced_pterostigmata()