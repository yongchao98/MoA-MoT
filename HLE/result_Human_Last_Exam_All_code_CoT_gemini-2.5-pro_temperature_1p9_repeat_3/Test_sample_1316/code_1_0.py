def explain_bridges_experiment():
    """
    Explains the non-disjunction event leading to an exceptional
    X0 red-eyed male in Drosophila.
    """
    # 1. Define Parental Genotypes
    female_parent = "XwXw"  # White-eyed female
    male_parent = "X+Y"     # Red-eyed (wild-type) male

    print("Analyzing the genetic cross to explain the unexpected offspring.")
    print(f"Parental Generation: White-eyed Female ({female_parent}) x Red-eyed Male ({male_parent})")
    print("-" * 60)

    # 2. Analyze the exceptional offspring
    exceptional_offspring = "X+0"
    phenotype = "Red-eyed Male"
    print(f"Unexpected Offspring Found: Genotype {exceptional_offspring}, Phenotype: {phenotype}")
    print("\nThis male must have received its X+ chromosome from the father, as the mother only has Xw alleles.")
    print("This means the mother must have contributed a gamete with no sex chromosome ('0').")
    print("-" * 60)

    # 3. Simulate the non-disjunction event
    event = "Non-disjunction of X chromosomes in female Meiosis I"
    print(f"Hypothesis: The causal event is '{event}'.")
    print("In this event, the female's homologous Xw chromosomes fail to separate,")
    print("producing two types of abnormal eggs: 'XwXw' and '0'.")
    print("-" * 60)

    # 4. Demonstrate the fertilization
    # Normal male gametes
    male_sperm_1 = "X+"
    male_sperm_2 = "Y"
    # Abnormal female gamete
    abnormal_egg = "0"

    print("The final fertilization event can be represented as an equation:")
    print(f"Paternal Gamete ({male_sperm_1}) + Maternal Gamete ({abnormal_egg})  =>  Zygote ({exceptional_offspring})")
    print("\nThis result perfectly matches the genotype of the unexpected red-eyed male.")
    print("Therefore, the non-disjunction occurred in female meiosis I.")


if __name__ == "__main__":
    explain_bridges_experiment()
