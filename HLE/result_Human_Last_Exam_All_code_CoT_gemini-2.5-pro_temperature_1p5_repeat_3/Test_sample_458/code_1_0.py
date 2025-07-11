def find_species_with_reduced_pterostigmata():
    """
    Identifies dragonfly species from a predefined list that are expected
    to have reduced pterostigmata based on their ecology.

    The primary factor is the trade-off between the pterostigma's role in
    preventing wing flutter during active flight and the aerodynamic drag
    it creates. Species that are long-distance migrants and rely heavily on
    gliding benefit more from reduced drag than from flutter suppression.

    The analysis identifies the following species as strong, gliding migrants:
    - 3: Macrodiplax balteata (Marl Pennant)
    - 4: Pantala flavescens (Globe Skimmer)
    - 8: Sympetrum corruptum (Variegated Meadowhawk)
    - 10: Tholymis tillarga (Coral-tailed Cloudwing)
    """

    # Indices of species known for long-distance gliding migration
    species_indices = [3, 4, 8, 10]

    # Format the output as a comma-separated string
    # The prompt requires outputting each number, which this loop does before the final result.
    print("Based on an analysis of their migratory and gliding ecologies, the following species indices correspond to those expected to have reduced pterostigmata:")
    for index in species_indices:
        print(f"Species Index: {index}")

    print("\nFinal comma-separated list:")
    # Using str.join to create the final comma-separated string from the integer list
    final_answer = ",".join(map(str, species_indices))
    print(final_answer)


if __name__ == "__main__":
    find_species_with_reduced_pterostigmata()