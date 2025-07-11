def explain_bridges_experiment():
    """
    Explains the genetic event leading to an X0 male Drosophila
    in the context of Bridges' experiments.
    """

    print("--- Analysis of an Exceptional Offspring in Bridges' Experiments ---")

    # 1. Define the observed offspring
    offspring_genotype = "X0"
    offspring_phenotype = "Male with red eyes and miniature wings"
    offspring_x_chromosome = "X(R,m)"  # (R) for red eyes, (m) for miniature wings

    print(f"Observed Offspring Genotype: {offspring_genotype}")
    print(f"Observed Offspring Phenotype: {offspring_phenotype}")
    print(f"This means its single X chromosome is: {offspring_x_chromosome}\n")

    # 2. Apply the context of Bridges' findings (female non-disjunction)
    print("--- Deducing the Causal Event ---")
    print("Context: Bridges' experiments demonstrated that exceptional offspring arise from non-disjunction in the female parent.")
    print("Conclusion: The '0' (nullo) gamete must be the egg from the mother, and the X chromosome must come from the father.\n")

    # 3. Determine the gametes involved
    paternal_gamete = offspring_x_chromosome
    maternal_gamete = "O (nullo-X egg)"

    print(f"Paternal Contribution: A normal sperm carrying the {paternal_gamete} chromosome.")
    print(f"Maternal Contribution: An abnormal egg with no X chromosome, which is a '{maternal_gamete}'.\n")

    # 4. Identify the meiotic error
    print("--- The Meiotic Error in the Female ---")
    print("A nullo-X egg is produced when the X chromosomes fail to separate during female meiosis.")
    print("This event is called non-disjunction.")
    print("Non-disjunction in Meiosis I is the failure of homologous X chromosomes to separate.")
    print("This process results in the formation of nullo-X eggs, along with XX eggs.\n")

    # 5. Reconstruct the fertilization event
    print("--- Final Equation of Fertilization ---")
    print(f"Maternal Gamete ({maternal_gamete}) + Paternal Gamete ({paternal_gamete}) --> Offspring ({offspring_genotype})")
    print("\nThis specific event, non-disjunction of the X chromosome in female meiosis I, explains the presence of the exceptional X0 male.")


explain_bridges_experiment()