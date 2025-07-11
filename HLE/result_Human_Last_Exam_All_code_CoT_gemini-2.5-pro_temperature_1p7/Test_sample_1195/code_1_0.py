def solve_genetics_problem():
    """
    Calculates the F2 phenotypic ratio for the given Drosophila cross.
    """
    print("### Step-by-Step Genetic Analysis ###\n")

    # Step 1: Define Parental (P) and F1 Generation
    print("Step 1: Determine the F1 generation from the parental cross.")
    print("Parental (P) Female Genotype: X(v)X(v); su-v/su-v")
    print("Parental (P) Male Genotype:   X(v)Y; su-v+/su-v+")
    print("All gametes from the female are X(v); su-v.")
    print("Gametes from the male are X(v); su-v+ or Y; su-v+.")
    print("\nResulting F1 Generation:")
    print("F1 Female Genotype: X(v)X(v); su-v+/su-v")
    print("F1 Male Genotype:   X(v)Y; su-v+/su-v")
    print("In the F1 generation, the suppressor gene is heterozygous (su-v+/su-v), so suppression does not occur. All F1 flies have vermilion eyes.\n")

    # Step 2: Analyze the F2 generation from the F1 cross
    print("Step 2: Analyze the F2 generation from the F1 cross (F1 x F1).")
    print("The cross is: X(v)X(v); su-v+/su-v  x  X(v)Y; su-v+/su-v\n")

    print("Analysis of the X-linked gene (X(v)X(v) x X(v)Y):")
    print("All parents in the F1 cross carry the X(v) allele.")
    print("Therefore, 100% of the F2 offspring will have the genotype for vermilion eyes (either X(v)X(v) or X(v)Y).")
    print("This means the final phenotype depends entirely on the autosomal suppressor gene.\n")

    # Step 3: Analyze the autosomal suppressor gene and calculate the final ratio
    print("Step 3: Calculate the phenotypic ratio based on the suppressor gene cross (su-v+/su-v x su-v+/su-v).")
    # In a standard monohybrid cross, the genotypic ratio is 1:2:1.
    # Total parts in the ratio = 1 + 2 + 1 = 4
    no_suppression_parts = 1 + 2  # for su-v+/su-v+ and su-v+/su-v
    suppression_parts = 1        # for su-v/su-v
    total_parts = 4

    print(f"The cross results in a genotypic ratio of 1 (su-v+/su-v+) : 2 (su-v+/su-v) : 1 (su-v/su-v).")
    print(f"Offspring with genotypes su-v+/su-v+ or su-v+/su-v will have NO suppression, resulting in vermilion eyes.")
    print(f"Offspring with genotype su-v/su-v will have suppression, resulting in wild-type eyes.\n")

    print("### Final F2 Phenotypic Ratio Calculation ###")
    print(f"Fraction with vermilion eyes (no suppression) = ({1} + {2}) / {total_parts} = {no_suppression_parts}/{total_parts}")
    print(f"Fraction with wild-type eyes (suppression) = {suppression_parts} / {total_parts}")

    print(f"\nThe final expected phenotypic ratio is {no_suppression_parts} vermilion : {suppression_parts} wild-type.")

# Run the function to display the results
solve_genetics_problem()
<<<B>>>