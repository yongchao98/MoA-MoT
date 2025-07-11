def find_species_with_reduced_pterostigmata():
    """
    This script identifies dragonfly species expected to have reduced pterostigmata
    based on their ecology.

    The principle is that long-distance migrants that utilize a gliding flight style
    (often called 'gliders') have a reduced need for the pterostigma, which acts
    as an anti-flutter device during active, high-speed flapping.

    The identified species are known long-distance gliders:
    - 3: Macrodiplax balteata (Marl Pennant)
    - 4: Pantala flavescens (Globe Skimmer)
    - 8: Sympetrum corruptum (Variegated Meadowhawk)
    - 10: Tholymis tillarga (Coral-tailed Cloudwing)
    """

    # The indices of the species expected to have reduced pterostigmata
    species_indices = [3, 4, 8, 10]

    # Print each number to form the final result string
    output = ",".join(map(str, species_indices))
    print(output)

find_species_with_reduced_pterostigmata()