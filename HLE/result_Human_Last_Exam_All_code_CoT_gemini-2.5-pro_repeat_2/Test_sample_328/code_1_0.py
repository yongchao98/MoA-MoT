def estimate_mating_age():
    """
    This function analyzes the life cycles of Tridactylophagus tartari and
    Periplaneta americana to estimate the average age of a mating male
    since eclosion.
    """

    # --- Species 1: Tridactylophagus tartari (Odd Beetle) ---
    species1_name = "Tridactylophagus tartari"
    age1_value = 3
    age1_unit = "days"
    explanation1 = (
        f"The adult {species1_name} has a very short lifespan (approx. 2 weeks) and does not feed. "
        "Its biological imperative is to mate as quickly as possible. "
        "Therefore, the best average estimate of a mating male's age is very short."
    )

    # --- Species 2: Periplaneta americana (American Cockroach) ---
    species2_name = "Periplaneta americana"
    age2_value = 6
    age2_unit = "months"
    explanation2 = (
        f"The adult {species2_name} lives for a year or more and can mate throughout this period. "
        "While it becomes sexually mature within a couple of weeks, a randomly observed mating "
        "male would, on average, be in the middle of its adult lifespan."
    )

    # --- Final Conclusion ---
    print("Ecological Observation Analysis:")
    print("-" * 50)

    print(f"Species: {species1_name}")
    print(explanation1)
    print(f"Estimated Age: {age1_value} {age1_unit}\n")

    print(f"Species: {species2_name}")
    print(explanation2)
    print(f"Estimated Age: {age2_value} {age2_unit}\n")

    print("-" * 50)
    print("Final combined estimate:")
    print(f"The best average estimate for {species1_name} is {age1_value} {age1_unit}, and for {species2_name} is {age2_value} {age2_unit}.")

# Execute the function to print the analysis
estimate_mating_age()
<<<C>>>