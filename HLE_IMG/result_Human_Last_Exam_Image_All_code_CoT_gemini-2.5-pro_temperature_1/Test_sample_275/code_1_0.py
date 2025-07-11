def calculate_splice_species():
    """
    Calculates the total number of splice species based on the provided RT-PCR results.

    The logic is to sum the number of unique species identified by each mutually exclusive
    forward primer class (E1F1, E1F2, E1F3).
    """

    # Number of species detected using the E1F1 forward primer.
    # Panels B, C, E, H, and I use this primer. Panel I is the most sensitive,
    # showing 4 distinct bands (450, 419, 378, and 187 bp).
    e1f1_species = 4

    # Number of species detected using the E1F2 forward primer (spanning 929^3465).
    # Panels D and F use this primer, both showing 1 band.
    e1f2_species = 1

    # Number of species detected using the E1F3 forward primer (spanning 929^3506).
    # Panel G uses this primer, showing 2 distinct bands (1,003 and 931 bp).
    e1f3_species = 2

    # The total number of species is the sum from each exclusive class.
    total_species = e1f1_species + e1f2_species + e1f3_species

    print("To find the total number of splice species, we sum the unique species identified by each primer class:")
    print(f"Number of species in E1F1 class (from Panel I): {e1f1_species}")
    print(f"Number of species in E1F2 class (from Panels D, F): {e1f2_species}")
    print(f"Number of species in E1F3 class (from Panel G): {e1f3_species}")
    print("\nFinal calculation:")
    # The final equation as requested.
    print(f"{e1f1_species} + {e1f2_species} + {e1f3_species} = {total_species}")
    print(f"\nTotal number of splice species found: {total_species}")

calculate_splice_species()
<<<7>>>