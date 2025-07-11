def explain_bridges_nondisjunction():
    """
    Explains the genetic cross and non-disjunction event from Bridges' experiments
    that leads to X0 male Drosophila.
    """

    # Let w = white-eye allele (recessive) and w+ = red-eye allele (dominant).
    # The gene for eye color is on the X chromosome.
    female_genotype = "X(w)/X(w)"  # White-eyed female
    male_genotype = "X(w+)/Y"     # Red-eyed male

    print("--- Calvin Bridges' Classic Cross ---")
    print(f"Parental (P) generation crosses:")
    print(f"  Female Parent Genotype: {female_genotype} (phenotype: white eyes)")
    print(f"  Male Parent Genotype:   {male_genotype} (phenotype: red eyes)")
    print("-" * 45)

    # --- Scenario 1: Normal Meiosis ---
    print("--- Outcome with Normal Meiosis ---")
    female_gamete_normal = "X(w)"
    male_gametes_normal = ["X(w+)", "Y"]
    print(f"The female produces one type of egg: {female_gamete_normal}")
    print(f"The male produces two types of sperm: {male_gametes_normal[0]} and {male_gametes_normal[1]}")

    # Expected offspring genotypes (F1)
    offspring1_equation = f"{female_gamete_normal} (egg) + {male_gametes_normal[0]} (sperm) = X(w+)/X(w)"
    offspring2_equation = f"{female_gamete_normal} (egg) + {male_gametes_normal[1]} (sperm) = X(w)/Y"

    print("\nExpected F1 Offspring:")
    print(f"  1. {offspring1_equation} -> Female, Red eyes")
    print(f"  2. {offspring2_equation} -> Male, White eyes")
    print("-" * 45)

    # --- Scenario 2: Non-disjunction in Female Meiosis I ---
    print("--- Outcome with Female Meiosis I Non-disjunction ---")
    print("In this rare event, the homologous X chromosomes in the female fail to separate.")
    
    # This leads to two types of abnormal eggs
    female_gametes_abnormal = ["X(w)X(w)", "0"] # '0' represents a gamete with no sex chromosome
    
    print(f"\nThis produces two types of abnormal eggs:")
    print(f"  1. Egg with both X chromosomes: {female_gametes_abnormal[0]}")
    print(f"  2. Egg with NO sex chromosome:  {female_gametes_abnormal[1]}")
    
    print("\nThe male produces normal sperm as before.")

    print("\nResulting 'Exceptional' Offspring:")
    # Exceptional Male: '0' egg fertilized by 'X(w+)' sperm from the father
    exceptional_male_equation = f"'{female_gametes_abnormal[1]}' (egg) + '{male_gametes_normal[0]}' (sperm) = X(w+)/0"
    
    print(f"  - Fertilization Equation: {exceptional_male_equation}")
    print("  - Resulting Phenotype: Male, Red eyes. This individual is an X0 male.")
    print("\nThis result precisely explains the 'unexpected' red-eyed male offspring.")
    print("The specific cause is the failure of the two X chromosomes to separate in the mother during meiosis.")
    print("This event is called Non-disjunction of the X chromosome in female meiosis I.")

explain_bridges_nondisjunction()