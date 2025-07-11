def explain_nondisjunction_in_drosophila():
    """
    This script explains the genetic events leading to an X0 male in Drosophila,
    as observed in Bridges' experiments.
    """
    
    # Define the parental genotypes and the observed unexpected offspring
    # In the classic cross, the male provides the red-eye allele (X+).
    # The X0 male's phenotype is determined by the father's X chromosome.
    father_contribution = "X-chromosome (with red-eye allele)"
    mother_contribution = "gamete with no sex chromosome ('0')"
    offspring_genotype = "X0"

    print("Step-by-step analysis of the unexpected X0 male offspring:")
    print("----------------------------------------------------------")
    
    print(f"1. The offspring's genotype is {offspring_genotype}. This means it has one X chromosome and no second sex chromosome.")
    
    print(f"2. To form this genotype, the zygote received an '{father_contribution}' from the father.")
    
    print(f"3. The zygote must have received a '{mother_contribution}' from the mother.")
    
    print("\n4. How can a female (XX) produce a '0' gamete (egg)?")
    print("   This occurs due to a meiotic error called non-disjunction, where chromosomes fail to separate.")
    
    print("\n5. Let's model the specific event: Non-disjunction in female Meiosis I.")
    print("   - In Meiosis I, the two homologous X chromosomes in the female fail to separate.")
    print("   - This results in the production of two types of eggs: eggs with two X's (XX) and eggs with zero X's (0).")
    
    print("\n6. The final fertilization equation is:")
    sperm_gamete = "X"
    egg_gamete = "0"
    zygote = "X0"
    # The instruction asked to "output each number in the final equation"
    # We will represent the gametes and zygote in the print statement.
    print(f"   Father's Sperm ({sperm_gamete}) + Mother's Egg ({egg_gamete}) --> Zygote ({zygote})")

    print("\nConclusion: The presence of X0 offspring with a paternal X chromosome directly indicates that non-disjunction of the X chromosomes occurred in the mother during egg formation. Non-disjunction in Meiosis I is the specific event that explains this.")

explain_nondisjunction_in_drosophila()