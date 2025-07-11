def estimate_mating_age():
    """
    This function provides the best average estimate for the age since eclosion
    for a male Tridactylophagus tartari and a male Periplaneta americana to be observed mating.
    """

    # --- Analysis ---
    # The term 'eclosion' is interpreted as hatching from the egg.
    
    # Species 1: Tridactylophagus tartari (a mite)
    # Mites in this group have extremely rapid life cycles. The time from hatching to
    # reaching sexual maturity can be as short as a few days.
    age_tridactylophagus = "three days"
    
    # Species 2: Periplaneta americana (American cockroach)
    # This species has a long development period. After eclosion (hatching),
    # it undergoes a lengthy nymphal stage before becoming a mature adult.
    # This nymphal period averages about six months.
    age_periplaneta = "six months"

    # --- Output ---
    print("Based on the typical life cycles of the two species:")
    
    # Final equation/statement components
    species_1_name = "Tridactylophagus tartari"
    species_2_name = "Periplaneta americana"
    
    print(f"The estimated age for a male {species_1_name} is: {age_tridactylophagus}")
    print(f"The estimated age for a male {species_2_name} is: {age_periplaneta}")

# Execute the function to print the result
estimate_mating_age()