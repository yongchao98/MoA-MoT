def bridges_experiment_simulation():
    """
    Simulates the chromosomal non-disjunction event described by Bridges
    that leads to an exceptional X0 male.
    """
    # Parent genotypes, based on the logic that the father provides the X chromosome.
    # The mother's genotype doesn't matter for the phenotype of this specific male,
    # but let's assume a common cross (e.g., wild-type female).
    mother_genotype = "X+/X+"
    father_genotype = "X^+m/Y"

    print(f"Parental Cross:")
    print(f"  Mother (Wild-type): {mother_genotype}")
    print(f"  Father (Red-eyed, Miniature wings): {father_genotype}\n")

    # --- Gamete Formation ---
    print("Gamete Formation:")

    # 1. Female Meiosis with Non-disjunction in Meiosis I
    # Homologous X chromosomes fail to separate, producing two types of eggs.
    female_nondisjunction_gamete_1 = "X+/X+" # Diploid egg
    female_nondisjunction_gamete_2 = "0"     # Nullo-X egg (no sex chromosome)
    print("  Mother undergoes Non-disjunction in Meiosis I, producing:")
    print(f"    - Diploid eggs: {female_nondisjunction_gamete_1}")
    print(f"    - Nullo-X eggs: {female_nondisjunction_gamete_2} <--- This causes the exception")


    # 2. Normal Male Meiosis
    male_gamete_1 = "X^+m"
    male_gamete_2 = "Y"
    print("\n  Father undergoes normal meiosis, producing:")
    print(f"    - Sperm 1: {male_gamete_1}")
    print(f"    - Sperm 2: {male_gamete_2}\n")

    # --- Fertilization producing the exceptional male ---
    print("Exceptional Fertilization Event:")
    egg = female_nondisjunction_gamete_2
    sperm = male_gamete_1
    offspring_genotype = sperm + egg # Order doesn't matter, but notation is usually X0

    print(f"  A nullo-X egg from the mother ({egg}) is fertilized by an X-bearing sperm from the father ({sperm}).\n")
    print("Resulting Offspring:")
    print(f"  Genotype: {offspring_genotype}")
    print("  Phenotype: Male (due to X0), Red eyes (from '+'), Miniature wings (from 'm')")

    print("\nConclusion: The presence of this X0 male indicates a non-disjunction of the X chromosome during female meiosis.")

bridges_experiment_simulation()