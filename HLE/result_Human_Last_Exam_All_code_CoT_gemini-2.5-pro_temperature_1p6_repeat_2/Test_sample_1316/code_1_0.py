def analyze_bridges_cross():
    """
    Simulates the genetic outcomes of Bridges' cross, focusing on the
    non-disjunction event that explains the exceptional X0 male.
    """
    # --- Parental Genotypes ---
    female_genotype = "XwXw"  # White-eyed female
    male_genotype = "X+Y"     # Red-eyed male

    print(f"Parental Cross: Female({female_genotype}) x Male({male_genotype})")
    print("-" * 50)
    print("Observed exceptional offspring: Red-eyed Male with X0 chromosomes.")
    print("This means his genotype is X+0.\n")

    print("To produce an X+0 individual:")
    print("1. The sperm from the father must have been: X+")
    print("2. The egg from the mother must have had no sex chromosome: 0 (nullo-X egg)")
    print("\nA nullo-X egg is created by a non-disjunction event in the female.\n")

    print("Let's analyze the outcome of 'Non-disjunction in Female Meiosis I' (Answer Choice A):")
    print("In this event, the two Xw chromosomes fail to separate.")
    # --- Gamete Formation ---
    # Normal male gametes
    male_gametes = ["X+", "Y"]
    # Abnormal female gametes due to Meiosis I non-disjunction
    female_abnormal_gametes = ["XwXw", "0"]

    print(f"This produces two types of abnormal eggs: {female_abnormal_gametes[0]} and {female_abnormal_gametes[1]}")
    print(f"The father produces two types of normal sperm: {male_gametes[0]} and {male_gametes[1]}\n")

    print("--- Fertilization 'Equations' and Resulting Offspring ---")

    # Combine abnormal eggs with normal sperm
    egg1 = female_abnormal_gametes[0] # "XwXw"
    egg2 = female_abnormal_gametes[1] # "0"
    sperm1 = male_gametes[0] # "X+"
    sperm2 = male_gametes[1] # "Y"

    # Equation 1: XwXw egg + X+ sperm
    zygote1 = egg1 + sperm1
    phenotype1 = "Dies (metafemale)"
    print(f"Egg({egg1}) + Sperm({sperm1}) -> Zygote({zygote1}). Phenotype: {phenotype1}")

    # Equation 2: XwXw egg + Y sperm
    zygote2 = egg1 + sperm2
    phenotype2 = "White-eyed Female (exceptional)"
    print(f"Egg({egg1}) + Sperm({sperm2}) -> Zygote({zygote2}). Phenotype: {phenotype2}")

    # Equation 3: 0 egg + X+ sperm
    zygote3 = sperm1 + egg2 # Swapped for conventional notation X0
    phenotype3 = "Red-eyed Male (exceptional, sterile)"
    print(f"Egg({egg2}) + Sperm({sperm1}) -> Zygote({zygote3}). Phenotype: {phenotype3}  <-- This matches the observation!")

    # Equation 4: 0 egg + Y sperm
    zygote4 = sperm2 + egg2
    phenotype4 = "Dies (non-viable)"
    print(f"Egg({egg2}) + Sperm({sperm2}) -> Zygote({zygote4}). Phenotype: {phenotype4}")

    print("-" * 50)
    print("Conclusion: The formation of the Red-eyed X0 Male is perfectly explained by")
    print("the fertilization of a nullo-X ('0') egg by an X+ sperm.")
    print("This specific type of egg is a direct result of non-disjunction of the X chromosome in female meiosis.")

# Run the analysis
analyze_bridges_cross()