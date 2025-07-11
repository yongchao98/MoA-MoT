def solve_bridges_experiment():
    """
    Analyzes the genetic scenario from Bridges' experiments to determine the cause of X0 males.
    
    The key steps are:
    1.  Identify the genotype of the exceptional offspring: X0 male.
    2.  Determine the gametes required to form this zygote: An X-bearing sperm and an egg with no sex chromosome ('0' egg).
    3.  In the classic cross (white-eyed female x red-eyed male), the exceptional male is red-eyed, meaning he got his X from his father.
        - Father (XY) must provide the X sperm.
        - Mother (XX) must provide the '0' egg.
    4.  A '0' egg is produced when the X chromosomes fail to separate during female meiosis (oogenesis).
    5.  This failure is called non-disjunction.
    6.  Non-disjunction of homologous chromosomes in Meiosis I in the female is the specific event that produces both '0' eggs (leading to X0 males) and 'XX' eggs (leading to the other exceptional offspring, XXY females).
    
    Therefore, the observation points to non-disjunction in female meiosis.
    """
    
    parent_female = "XX"
    parent_male = "XY"
    
    exceptional_offspring = "X0"
    
    # To get X0, we need one gamete with 'X' and one with '0' (no sex chromosome)
    # The X must come from the father in the classic experiment to explain the phenotype.
    sperm_gamete = "X"
    egg_gamete = "0"
    
    print(f"To produce an exceptional {exceptional_offspring} male, the zygote must be formed from a sperm and an egg.")
    print(f"Based on the inheritance pattern of X-linked traits in Bridges' experiment, the sperm must contribute the 'X' chromosome.")
    print(f"Therefore, the egg from the female parent ({parent_female}) must have contributed no sex chromosome, being a '{egg_gamete}' gamete.")
    
    cause = "Non-disjunction of the X chromosome in the female parent during meiosis"
    
    print(f"An egg with no sex chromosome is the result of: {cause}.")
    print("This can happen in Meiosis I (failure of homologous chromosomes to separate) or Meiosis II (failure of sister chromatids to separate).")
    print("Non-disjunction in Meiosis I is the primary event that also explains the concurrent appearance of exceptional XXY females.")
    print("\nConclusion: The presence of X0 males indicates a non-disjunction event of the X chromosome during female gametogenesis.")

# Execute the function to print the explanation.
solve_bridges_experiment()