import fractions

def solve_genetics_cross():
    """
    Solves the Drosophila eye color genetics problem by calculating the F2 phenotypic ratio.
    """

    print("Step 1: Define Parental (P) Generation Genotypes")
    # v = vermilion allele (X-linked, recessive)
    # + = wild-type allele for suppressor gene (autosomal, dominant)
    # su-v = suppressor of vermilion allele (autosomal, recessive)
    # An X^v fly that is also su-v/su-v has its vermilion color suppressed, resulting in a wild-type phenotype.
    p_female_genotype = "X(v)X(v); su-v/su-v"
    p_male_genotype = "X(v)Y; +/+"
    print(f"Parental Female: {p_female_genotype}")
    print(f"Parental Male:   {p_male_genotype}\n")

    print("Step 2: Determine F1 Generation Genotypes")
    # Female gametes are all X(v); su-v
    # Male gametes are X(v); + and Y; +
    f1_female_genotype = "X(v)X(v); +/su-v" # From (X(v); su-v) x (X(v); +)
    f1_male_genotype = "X(v)Y; +/su-v"   # From (X(v); su-v) x (Y; +)
    print("All F1 females will be heterozygous for the suppressor gene and have vermilion eyes.")
    print(f"F1 Female Genotype: {f1_female_genotype} (Phenotype: Vermilion)")
    print("All F1 males will be heterozygous for the suppressor gene and have vermilion eyes.")
    print(f"F1 Male Genotype:   {f1_male_genotype} (Phenotype: Vermilion)\n")

    print("Step 3: Determine F2 Generation from F1 Cross")
    print(f"F2 Cross: {f1_female_genotype} (Female) x {f1_male_genotype} (Male)")
    
    # Analysis of the autosomal suppressor gene cross: +/su-v x +/su-v
    # This standard monohybrid cross gives a genotypic ratio of 1:2:1
    prob_plus_plus = 0.25      # Genotype: +/+
    prob_plus_su_v = 0.50      # Genotype: +/su-v
    prob_su_v_su_v = 0.25      # Genotype: su-v/su-v
    
    print("\nSince all F1 parents have the X(v) gene, all F2 offspring will also have the genetic potential for vermilion eyes.")
    print("Therefore, the final phenotype depends solely on the autosomal suppressor gene (`su-v`).\n")

    print("Step 4: Calculate the F2 Phenotypic Ratio")
    
    # Phenotypes are determined by whether the suppression is active.
    # Suppression is INACTIVE if at least one '+' allele is present.
    # Suppression is ACTIVE if the genotype is homozygous recessive 'su-v/su-v'.
    
    prob_vermilion = prob_plus_plus + prob_plus_su_v
    prob_wild_type = prob_su_v_su_v
    
    print("- Flies with genotypes +/+ or +/su-v will NOT have their vermilion color suppressed.")
    print(f"  Fraction of vermilion flies = {fractions.Fraction(prob_plus_plus).limit_denominator()} (+/+) + {fractions.Fraction(prob_plus_su_v).limit_denominator()} (+/su-v) = {fractions.Fraction(prob_vermilion).limit_denominator()}")

    print("- Flies with genotype su-v/su-v WILL have their vermilion color suppressed, resulting in wild-type eyes.")
    print(f"  Fraction of wild-type flies = {fractions.Fraction(prob_wild_type).limit_denominator()} (su-v/su-v)")

    # The ratio of numerators
    numerator_wild_type = fractions.Fraction(prob_wild_type).limit_denominator().numerator
    numerator_vermilion = fractions.Fraction(prob_vermilion).limit_denominator().numerator

    print("\n--- Final F2 Phenotypic Ratio ---")
    print(f"The final predicted ratio is:")
    print(f"{numerator_vermilion} vermilion : {numerator_wild_type} wild-type")
    print("Expressed as fractions of the total F2 population:")
    print(f"Eye Color Phenotype Breakdown:")
    print(f"  Wild-type: {fractions.Fraction(prob_wild_type).limit_denominator()}")
    print(f"  Vermilion: {fractions.Fraction(prob_vermilion).limit_denominator()}")

solve_genetics_cross()
<<<B>>>