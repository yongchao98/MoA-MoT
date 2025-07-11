def analyze_bridges_experiment():
    """
    This script analyzes the genetic cross from Bridges' experiment to pinpoint
    the specific non-disjunction event responsible for the exceptional offspring.
    """

    # Step 1: Define the parental cross and the exceptional offspring in question.
    # The mother is white-eyed, so her genotype is X(w)X(w).
    # The father is red-eyed, so his genotype is X(R)Y.
    # The unexpected offspring is a red-eyed male with an X0 chromosomal makeup.
    mother_genotype = "X(w)X(w)"
    father_genotype = "X(R)Y"
    exceptional_male_genotype = "X(R)0"

    print("--- Analysis of Bridges' Experiment ---")
    print(f"Parental Cross: White-eyed Female ({mother_genotype}) x Red-eyed Male ({father_genotype})")
    print(f"Exceptional Offspring: Red-eyed Male with Genotype {exceptional_male_genotype}")
    print("-" * 35 + "\n")

    # Step 2: Deduce the gametes that formed this exceptional offspring.
    # The male's red eyes (R) must come from the father's X(R) chromosome.
    # The male's X0 makeup means he received no sex chromosome from his mother.
    father_gamete = "X(R)"
    mother_gamete = "0" # '0' represents a gamete with no sex chromosome (a nullo-X egg).

    print("--- Deducing the Gametes ---")
    print(f"1. The red-eye allele (R) must come from the father ({father_genotype}). So, the father contributed an '{father_gamete}' sperm.")
    print(f"2. The X0 genotype means the mother ({mother_genotype}) contributed a gamete with no X chromosome, which we represent as '{mother_gamete}'.")
    print("-" * 35 + "\n")

    # Step 3: Reconstruct the fertilization event as an equation.
    print("--- The Final Equation for the Exceptional Male ---")
    # The prompt asks to output each part of the final equation.
    print(f"Father's Gamete: {father_gamete}")
    print(f"Mother's Gamete: {mother_gamete}")
    print(f"Resulting Zygote: {exceptional_male_genotype}")
    print(f"Equation: {father_gamete} + {mother_gamete} = {exceptional_male_genotype}")
    print("-" * 35 + "\n")

    # Step 4: Determine the meiotic error that produced the mother's abnormal gamete.
    print("--- Identifying the Causal Event ---")
    print(f"The mother ({mother_genotype}) produced a '{mother_gamete}' egg. This requires a non-disjunction event during her meiosis.")
    print("\nLet's consider the options:")
    print("A. Non-disjunction in female meiosis I:")
    print("   - The two homologous X(w) chromosomes fail to separate.")
    print("   - This event produces two kinds of eggs: X(w)X(w) and '0'.")
    print("   - This correctly explains the origin of the '0' egg.")
    print("\nB. Non-disjunction in female meiosis II:")
    print("   - Sister chromatids of one X(w) chromosome fail to separate.")
    print("   - This also produces '0' eggs, but it is the failure of homologous chromosomes to separate (Meiosis I) that is the classic explanation for Bridges' results, as it provides a single mechanism for both the '0' eggs and the 'X(w)X(w)' eggs (which create the exceptional white-eyed females).")
    print("\nConclusion: The most direct and comprehensive explanation is non-disjunction during female meiosis I.")

analyze_bridges_experiment()