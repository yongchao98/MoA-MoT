def count_splice_species():
    """
    This function calculates the total number of unique splice species based on the provided RT-PCR data.
    """
    print("Analyzing the RT-PCR results to identify unique splice species:")

    # Panel B identifies two initial early splice variants.
    species_from_B = 2
    print(f"- Panel B (E1F1/E4R) reveals {species_from_B} distinct early splice species.")

    # Panel D uses a splice-junction-specific primer (E1F2), identifying a new species.
    species_from_D = 1
    print(f"- Panel D (E1F2/E5R) reveals {species_from_D} new early splice species, defined by the 929^3465 splice junction.")

    # Panel G uses another splice-junction-specific primer (E1F3), identifying two more new species.
    species_from_G = 2
    print(f"- Panel G (E1F3/L2R) reveals {species_from_G} new early splice species, defined by the 929^3506 splice junction.")

    # Panels H and I together identify all late splice variants.
    species_from_HI = 4
    print(f"- Panels H and I (E1F1/L1R) together reveal {species_from_HI} distinct late splice species.")
    
    # The total number of species is the sum of unique species found in each key experiment.
    total_species = species_from_B + species_from_D + species_from_G + species_from_HI

    print("\nCalculating the total number of unique splice species:")
    # The final equation prints each number explicitly as requested
    print(f"Total Species = {species_from_B} (from B) + {species_from_D} (from D) + {species_from_G} (from G) + {species_from_HI} (from H/I)")
    print(f"Total Species = {total_species}")

count_splice_species()