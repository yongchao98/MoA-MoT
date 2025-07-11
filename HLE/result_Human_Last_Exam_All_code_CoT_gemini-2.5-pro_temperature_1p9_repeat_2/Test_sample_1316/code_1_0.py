def analyze_bridges_nondisjunction():
    """
    This program explains the chromosomal event behind the findings in
    Calvin Bridges' experiments with Drosophila.
    """
    
    # 1. Define the classic parental cross setup.
    # The female has the white-eye allele (w) on both X's.
    # The male has the red-eye allele (+) on his X.
    # The question also mentions miniature wings, but eye color is sufficient to demonstrate the principle.
    female_genotype = "X(w)X(w)"
    male_genotype = "X(+)Y"
    
    print("--- Analysis of Bridges' Exceptional Offspring ---")
    print(f"Parental Cross: White-eyed Female ({female_genotype}) x Red-eyed Male ({male_genotype})\n")

    # 2. Describe the unexpected offspring in the question.
    print("Observed Anomaly: A male fly with red eyes and an X0 chromosomal makeup.")
    print("This means his genotype is X(+)0.\n")

    # 3. Determine the origin of the gametes that formed the X0 male.
    print("To form an X(+)0 individual, two gametes must fuse:")
    print(" - One gamete provides the X(+) chromosome.")
    print(" - One gamete provides NO sex chromosome (a '0' gamete).\n")
    print(f"The father ({male_genotype}) is the only parent with an X(+) chromosome.")
    print("Therefore, the X(+) came from the father's sperm.")
    print("This implies the '0' gamete (an egg with no sex chromosome) came from the mother.\n")
    
    # 4. Explain the mechanism: Non-disjunction in Female Meiosis I.
    print("--- The Causal Event: Non-disjunction in Female Meiosis I ---")
    print(f"The mother's genotype is {female_genotype}.")
    print("During Meiosis I, the homologous X chromosomes fail to separate.")
    print("This single error produces two kinds of abnormal eggs:")
    print(f" - An 'XX' egg containing both chromosomes: {female_genotype}")
    print(f" - A 'nullo-X' egg containing no sex chromosome: 0\n")

    # 5. Show the resulting fertilization "equations".
    print("--- Reconstructing the Exceptional Offspring ---")
    
    # The equation for the exceptional male
    abnormal_egg_0 = "0"
    normal_sperm_X_plus = "X(+)"
    print(f"Equation for the exceptional male:")
    print(f"   Abnormal Egg ({abnormal_egg_0}) + Normal Sperm ({normal_sperm_X_plus})  -->  {normal_sperm_X_plus}{abnormal_egg_0}")
    print("   Result: An X0 male with red eyes (paternal phenotype). This matches the observation.\n")
    
    # The equation for the exceptional female (the other half of Bridges' proof)
    abnormal_egg_XX = female_genotype
    normal_sperm_Y = "Y"
    print("Equation for the exceptional female (Bridges' other key finding):")
    print(f"   Abnormal Egg ({abnormal_egg_XX}) + Normal Sperm ({normal_sperm_Y})  -->  {abnormal_egg_XX}{normal_sperm_Y}")
    print("   Result: An XXY female with white eyes (maternal phenotype).\n")

    print("Conclusion: Non-disjunction of the X chromosome in female meiosis I is the single")
    print("event that accounts for the production of the nullo-X egg required to form the X0 male.")

analyze_bridges_nondisjunction()
<<<A>>>