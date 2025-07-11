def estimate_mating_age():
    """
    Analyzes the life cycles of Tridactylophagus tartari and Periplaneta americana
    to estimate the average age of a mating male.
    """
    
    # Species 1: Tridactylophagus tartari (a blister beetle)
    # This species has a short adult lifespan and mates early.
    beetle_species = "Tridactylophagus tartari"
    beetle_estimated_age = "three days"
    
    # Species 2: Periplaneta americana (American cockroach)
    # This species has a long adult lifespan and mates throughout its life.
    cockroach_species = "Periplaneta americana"
    cockroach_estimated_age = "six months"

    print("Estimating the average age of mating males for two different insect species.")
    print("-" * 60)
    
    print(f"\nAnalysis for {beetle_species}:")
    print("Blister beetles like this one typically have short adult lifespans (weeks to a few months).")
    print("They focus on reproduction shortly after emerging, often mating within the first few days.")
    print("Therefore, a 'best average estimate' for a mating male's age would be very young.")
    print(f"Conclusion for {beetle_species}: The estimated age is {beetle_estimated_age}.")
    
    print(f"\nAnalysis for {cockroach_species}:")
    print("The American cockroach has a long adult lifespan, often over a year.")
    print("While it reaches sexual maturity in about a week, it can mate throughout its life.")
    print("The 'average estimate' of a mating male's age would therefore be a significant fraction of its total lifespan, not its minimum mating age.")
    print(f"Conclusion for {cockroach_species}: The estimated age is {cockroach_estimated_age}.")
    
    print("-" * 60)
    print("\nFinal Answer:")
    print(f"The best average estimate for the {beetle_species} male is {beetle_estimated_age}.")
    print(f"The best average estimate for the {cockroach_species} male is {cockroach_estimated_age}.")
    
estimate_mating_age()
<<<C>>>