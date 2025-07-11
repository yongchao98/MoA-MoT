def bridges_experiment_simulation():
    """
    Simulates the genetic cross from Bridges' experiments on Drosophila
    to explain exceptional offspring.
    """

    # --- Step 1: Define Parental Genotypes ---
    female_parent = "XwXw (White-eyed female)"
    male_parent = "X+Y (Red-eyed male)"
    print("--- Bridges' Experiment Simulation ---")
    print(f"Parental Cross: {female_parent} x {male_parent}\n")

    # --- Step 2: Normal Meiosis and Offspring ---
    print("--- Part A: Normal Meiosis ---")
    print(f"The {female_parent} normally produces eggs with one genotype: Xw")
    print(f"The {male_parent} normally produces sperm with two genotypes: X+ and Y")
    print("\nExpected Offspring (Normal Inheritance):")
    print("1. Xw (egg) + X+ (sperm) = X+Xw (Red-eyed female)")
    print("2. Xw (egg) + Y (sperm)  = XwY (White-eyed male)")
    print("-" * 35)

    # --- Step 3 & 4: Simulate Meiosis I Non-disjunction and Analyze Outcome ---
    print("\n--- Part B: Non-Disjunction in Female Meiosis I ---")
    print("Hypothesis: Non-disjunction of X chromosomes occurs during Meiosis I in the female.")
    print("In this event, the two Xw chromosomes fail to separate.")
    print("This leads to two types of abnormal eggs being produced:")
    print("1. An egg with BOTH X chromosomes: XwXw")
    print("2. An egg with NO X chromosome: O (nullo-X)")

    print("\n--- Part C: Fertilization with Abnormal Eggs ---")
    print("These abnormal eggs are fertilized by the normal sperm from the male (X+ and Y):")

    # Fertilization Combination 1: XwXw egg
    print("\nFertilizing the 'XwXw' egg:")
    print("Equation 1: XwXw (from female) + Y (from male) = XwXwY (Exceptional white-eyed female)")
    print("Equation 2: XwXw (from female) + X+ (from male) = X+XwXw (Metafemale, usually dies)")

    # Fertilization Combination 2: O (nullo-X) egg
    print("\nFertilizing the 'O' (nullo-X) egg:")
    print("This combination produces the male offspring described in the question.")
    print("\nFinal Equation for the unexpected male:")
    # This is the key calculation the user requested.
    print(">>> O (from female) + X+ (from male) = X+O <<<")
    print("\nOffspring Analysis:")
    print("Genotype: X+O")
    print("Phenotype: The X+ from the father confers red eyes. The O from the mother results in a sterile male.")
    print("This perfectly matches the unexpected red-eyed male offspring with X0 chromosomes.")

    print("\nConclusion:")
    print("The presence of X+O males indicates that the mother produced 'O' (nullo-X) eggs.")
    print("This specific error is caused by the non-disjunction of X chromosomes in the female during meiosis.")
    print("This event, specifically in Meiosis I, also produces the reciprocal XwXw egg, which accounts for the exceptional white-eyed females (XwXwY) that Bridges also observed.")

# Run the simulation
bridges_experiment_simulation()