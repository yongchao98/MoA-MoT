def explain_bridges_experiment():
    """
    Explains the non-disjunction event in Bridges' Drosophila experiment
    that leads to exceptional red-eyed X0 males.
    """
    # Define parental genotypes and phenotypes
    # In Drosophila, the white-eye allele (w) is recessive to the red-eye allele (R),
    # and the gene is located on the X chromosome.
    female_genotype = "XwXw"
    female_phenotype = "white-eyed"
    male_genotype = "X(R)Y"
    male_phenotype = "red-eyed"

    print("--- Analysis of Bridges' Drosophila Experiment ---")
    print("\nStep 1: The Parental Cross")
    print(f"Parental Female: Genotype {female_genotype}, Phenotype: {female_phenotype}")
    print(f"Parental Male:   Genotype {male_genotype}, Phenotype: {male_phenotype}")

    print("\nStep 2: Expected Offspring from Normal Meiosis")
    print("A normal male (X(R)Y) produces X(R) sperm and Y sperm.")
    print("A normal female (XwXw) produces only Xw eggs.")
    print("Expected Female Offspring: X(R) sperm + Xw egg -> X(R)Xw (Red-eyed)")
    print("Expected Male Offspring:   Y sperm + Xw egg -> XwY (White-eyed)")

    print("\nStep 3: The Unexpected Observation")
    exceptional_male_genotype = "X(R)0"
    exceptional_male_phenotype = "red-eyed"
    print(f"An unexpected male was found with phenotype '{exceptional_male_phenotype}' and genotype '{exceptional_male_genotype}'.")
    print("(The '0' indicates the absence of a second sex chromosome).")

    print("\nStep 4: Deducing the Gametes for the Unexpected Male")
    print(f"To create an {exceptional_male_genotype} individual, we can trace the gametes:")
    father_gamete = "X(R) sperm"
    mother_gamete = "0 (nullo-X) egg"
    print(f"  - The red-eye allele (R) must come from the father, so he contributed an {father_gamete}.")
    print(f"  - To have no second sex chromosome, the mother must have contributed an egg with no X chromosome, an {mother_gamete}.")
    
    print("\nFinal Equation for the Unexpected Male:")
    print(f"Gamete from Father ({father_gamete}) + Gamete from Mother ({mother_gamete}) -> Offspring ({exceptional_male_genotype})")

    print("\nStep 5: Identifying the Causal Event")
    print(f"The mother ({female_genotype}) produced an abnormal '{mother_gamete}'.")
    print("This is caused by a failure of the X chromosomes to separate during her meiosis, an event called non-disjunction.")
    print("Specifically, when the two homologous X chromosomes fail to separate during Meiosis I, two types of abnormal eggs are produced: XwXw eggs and nullo-X (0) eggs.")
    print("This event, called primary non-disjunction, perfectly explains the formation of the nullo-X egg required.")

    print("\nConclusion:")
    print("The presence of an X(R)0 male proves that the mother produced a nullo-X egg, which resulted from the non-disjunction of her homologous X chromosomes during Meiosis I.")

# Execute the explanation
explain_bridges_experiment()