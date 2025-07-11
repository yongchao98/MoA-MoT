def solve_drosophila_genetics():
    """
    This function calculates the expected phenotypic ratio for the described
    Drosophila cross.
    """
    print("This script calculates the F2 phenotypic ratio for a dihybrid cross in Drosophila.")
    print("The cross involves an X-linked gene (vermilion, v) and an autosomal suppressor (su-v).\n")

    # Parental (P) Generation Genotypes
    p_female_genotype = "X(v)X(v); su-v/su-v"
    p_male_genotype = "X(v)Y; +/+"
    
    print("--- Step 1: Parental (P) Generation ---")
    print(f"Female Genotype: {p_female_genotype}")
    print(f"Male Genotype: {p_male_genotype}\n")
    
    print("--- Step 2: F1 Generation ---")
    # All F1 offspring are heterozygous for the suppressor gene (+/su-v) and have the vermilion allele (v).
    # Their phenotype is vermilion because the dominant '+' allele is present, preventing suppression.
    f1_genotype_info = "All F1 flies are heterozygous for the suppressor gene (+/su-v) and have the vermilion allele (v)."
    print(f1_genotype_info)
    print("F1 Phenotype: All Vermilion\n")

    print("--- Step 3: F2 Generation (from F1 x F1 cross) ---")
    # The F1 cross is: X(v)X(v); +/su-v  x  X(v)Y; +/su-v
    # For the X-linked gene: All F2 flies inherit an X(v) chromosome.
    # For the autosomal gene (+/su-v x +/su-v), the genotypic ratio is 1 (+/+): 2 (+/su-v): 1 (su-v/su-v).
    
    # Phenotype depends on the suppressor gene's genotype.
    # Vermilion phenotype: Occurs with +/+ or +/su-v genotypes.
    # Wild-type phenotype: Occurs with su-v/su-v genotype (suppression).

    vermilion_num = 3
    wildtype_num = 1
    denominator = 4
    
    print("In the F2 generation, the phenotype depends on the autosomal suppressor gene:")
    print("- Flies with genotypes +/+ or +/su-v are not suppressed and will be vermilion.")
    print("- Flies with genotype su-v/su-v are suppressed and will be wild-type.\n")
    
    print("--- Step 4: Final F2 Phenotypic Ratio ---")
    print(f"Fraction of vermilion flies = 1/4 (+/+) + 2/4 (+/su-v) = {vermilion_num}/{denominator}")
    print(f"Fraction of wild-type flies = 1/4 (su-v/su-v) = {wildtype_num}/{denominator}")
    print("\nThe final predicted ratio is:")
    print(f"{vermilion_num}/{denominator} vermilion : {wildtype_num}/{denominator} wild-type")

solve_drosophila_genetics()
<<<B>>>