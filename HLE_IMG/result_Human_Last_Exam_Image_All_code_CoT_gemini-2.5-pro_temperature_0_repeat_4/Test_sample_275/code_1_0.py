def count_splice_species():
    """
    Calculates the total number of splice species based on the provided RT-PCR data.
    """
    # Number of species from Panel B (E1F1/E4R)
    species_from_B = 2  # Bands at 190 bp and 159 bp

    # Number of new species from Panel D (E1F2/E5R)
    # This panel identifies a species with the 929^3465 splice.
    species_from_D = 1  # Band at 551 bp

    # Number of new species from Panel E (E1F1/L2R)
    # This panel identifies a distinct early transcript.
    species_from_E = 1 # Band at 1005 bp

    # Number of new species from Panel G (E1F3/L2R)
    # This panel identifies two species with the 929^3506 splice.
    species_from_G = 2  # Bands at 1003 bp and 931 bp

    # Number of new species from Panel I (E1F1/L1R)
    # This panel identifies four distinct late transcripts.
    species_from_I = 4  # Bands at 450, 419, 378, and 187 bp

    # Total number of species is the sum of unique species from each informative experiment.
    total_species = species_from_B + species_from_D + species_from_E + species_from_G + species_from_I

    print("To find the total number of splice species, we sum the unique variants identified in each informative experiment:")
    print(f"Panel B revealed {species_from_B} species.")
    print(f"Panel D revealed {species_from_D} new species.")
    print(f"Panel E revealed {species_from_E} new species.")
    print(f"Panel G revealed {species_from_G} new species.")
    print(f"Panel I revealed {species_from_I} new (late) species.")
    print("\nThe final calculation is:")
    print(f"{species_from_B} + {species_from_D} + {species_from_E} + {species_from_G} + {species_from_I} = {total_species}")
    print(f"\nTherefore, a total of {total_species} splice species can be found in this image.")

count_splice_species()