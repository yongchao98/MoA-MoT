def explain_resistance_acquisition():
    """
    This script explains how a bacterium with a stable genome can acquire
    drug resistance as rapidly as one with frequent lateral gene transfer.
    """

    # --- Scenario Setup ---
    # Let's assign numerical values to represent the 'fitness' of each bacterial type.
    # Fitness is a measure of reproductive success. A higher value means it grows faster.
    fitness_wild_type = 1.0
    fitness_resistant_with_cost = 0.8  # Resistance often comes with a biological cost
    fitness_compensated_resistant = 1.5 # A second mutation can compensate for the cost and even provide a new advantage

    print("Analyzing two models of drug resistance acquisition:")
    print("-" * 50)

    # --- Bacterium 1: Lateral Gene Transfer (LGT) ---
    print("Model 1: Bacterium with Lateral Gene Transfer")
    print("This bacterium can quickly acquire pre-existing resistance genes from other bacteria. This is a very fast method of spread.\n")

    # --- Bacterium 2: Mutation and Selection ---
    print("Model 2: Bacterium with a Stable Genome (Mutation-driven)")
    print("This bacterium must rely on random mutations. This is typically a slower process.")
    print("How can it keep pace with Model 1? Let's look at the steps:")

    # Step 1: The initial resistance mutation
    print("\nStep A: A rare mutation provides drug resistance.")
    print(f"However, this mutation often has a fitness cost, making the bacterium less healthy in a drug-free environment.")
    print(f"The equation of fitness change is: {fitness_wild_type} -> {fitness_resistant_with_cost}. The bacterium is now less competitive.")

    # Step 2: The compensatory mutation
    print("\nStep B: A second, 'compensatory' mutation occurs.")
    print(f"This second mutation not only fixes the problem from the first one but can dramatically increase overall fitness.")
    print(f"The equation of fitness change is: {fitness_resistant_with_cost} -> {fitness_compensated_resistant}. This strain is now highly competitive.")

    # Step 3: Rapid spread and cross-resistance
    print("\nStep C: The highly fit, resistant strain takes over.")
    print(f"Because its fitness ({fitness_compensated_resistant}) is much higher than the original wild-type ({fitness_wild_type}), it rapidly dominates the population through natural selection.")
    print("Furthermore, if this set of mutations also provides 'cross-resistance' (resistance to multiple drugs), it represents a huge leap in adaptation.")

    # --- Conclusion ---
    print("\n" + "-" * 50)
    print("Conclusion: The combination of a rare resistance mutation followed by a compensatory mutation that greatly increases fitness allows the mutation-driven evolution to be surprisingly fast. Adding cross-resistance further accelerates the perceived pace of adaptation. This makes choice B the most complete explanation.")


explain_resistance_acquisition()
<<<B>>>