def explain_bridges_experiment():
    """
    Explains the genetic cross and non-disjunction event that produces an X0 male in Drosophila.
    """
    print("Analyzing the genetic event in Misión Bridges' experiment:")
    print("-" * 60)

    # The offspring's phenotype is determined by the father's X chromosome.
    # Let's denote the allele for red eyes as '+' and miniature wings as 'm'.
    # The father's X-bearing sperm would be 'X+m'.
    # The mother's genotype is irrelevant for the non-disjunction mechanism itself.

    print("Step 1: Define the key offspring.")
    offspring_genotype = "X0"
    offspring_phenotype = "Male with red eyes and miniature wings"
    print(f"The unexpected offspring is an {offspring_genotype} male. In Drosophila, this individual is sterile.")

    print("\nStep 2: Determine the origin of the single X chromosome.")
    print("This male must inherit his single X chromosome from his father.")
    print("This means the mother must have produced an egg with NO sex chromosome.")
    maternal_gamete = "O" # 'O' represents the absence of a sex chromosome.
    paternal_gamete = "X+m" # Sperm carrying the father's X chromosome.

    print(f"\nStep 3: Identify the cause of the abnormal maternal gamete.")
    print(f"The production of a nullo-X ('{maternal_gamete}') egg is due to a failure of X chromosomes to separate during female meiosis.")
    print("This error is called non-disjunction.")
    print("This non-disjunction event is the direct cause of the unexpected offspring.")

    print("\nStep 4: Show the final fertilization equation.")
    print("The final equation for the formation of this unexpected male is:")
    # Using f-string to print each component of the equation clearly.
    print(f"Maternal Egg({maternal_gamete}) + Paternal Sperm({paternal_gamete}) ---> Zygote({offspring_genotype})")

    print("\nConclusion:")
    print("The event must be non-disjunction of the X chromosome in the female during meiosis.")
    print("Bridges' 'primary non-disjunction' refers to the failure of homologous chromosomes to separate, which occurs during Meiosis I.")
    print("Therefore, this points to non-disjunction of the X chromosome in female meiosis I.")

explain_bridges_experiment()