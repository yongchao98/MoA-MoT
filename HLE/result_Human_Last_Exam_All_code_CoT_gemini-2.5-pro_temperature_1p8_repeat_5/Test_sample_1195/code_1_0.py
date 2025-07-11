def solve_genetics_cross():
    """
    Calculates and explains the F2 phenotypic ratio for a Drosophila cross
    involving an X-linked gene and an autosomal suppressor gene.
    """
    # Step 1: Define Parental (P) generation
    p_female_genotype = "X(v)X(v); su-v/su-v"
    p_male_genotype = "X(v)Y; su-v+/su-v+"
    print("--- Step 1: Parental (P) Generation ---")
    print(f"Female Genotype: {p_female_genotype} (Phenotype: Wild-type due to suppressor)")
    print(f"Male Genotype:   {p_male_genotype} (Phenotype: Vermilion)\n")

    # Step 2: Determine F1 generation
    f1_female_genotype = "X(v)X(v); su-v+/su-v"
    f1_male_genotype = "X(v)Y; su-v+/su-v"
    print("--- Step 2: F1 Generation (from P cross) ---")
    print(f"All F1 Female Genotypes: {f1_female_genotype}")
    print(f"All F1 Male Genotypes:   {f1_male_genotype}")
    print("F1 Phenotype: All F1 flies are vermilion because the suppressor allele (su-v) is recessive.\n")

    # Step 3: Analyze the F1 intercross (F1 x F1) to get the F2 generation
    print("--- Step 3: F2 Generation Analysis (from F1 x F1 cross) ---")
    print("The F1 cross is: X(v)X(v); su-v+/su-v (female) x X(v)Y; su-v+/su-v (male).")
    print("Note: All F2 offspring receive a vermilion allele 'X(v)' from the F1 mother.")
    print("Therefore, the final eye color is determined solely by the 'su-v' gene locus.\n")
    
    # Step 4: Analyze the autosomal suppressor gene cross
    print("--- Step 4: Analyzing the 'su-v' Gene Cross ---")
    print("The cross is su-v+/su-v x su-v+/su-v. This gives a 3:1 phenotypic ratio.")
    
    # Define probabilities based on the monohybrid cross
    total_parts = 4
    # 3/4 of progeny will have genotype su-v+/_ (non-suppressing)
    non_suppressing_parts = 3  
    # 1/4 of progeny will have genotype su-v/su-v (suppressing)
    suppressing_parts = 1      

    # Step 5: Calculate and display the final F2 phenotypic ratio
    print("\n--- Step 5: Final F2 Phenotypic Ratio ---")
    print(f"Fraction with non-suppressing genotype (su-v+/_), resulting in vermilion eyes: {non_suppressing_parts}/{total_parts}")
    print(f"Fraction with suppressing genotype (su-v/su-v), resulting in wild-type eyes: {suppressing_parts}/{total_parts}")
    
    print("\nFinal Equation for Phenotypic Ratio [vermilion : wild-type]:")
    print(f"{non_suppressing_parts}/{total_parts} : {suppressing_parts}/{total_parts}")
    
    print(f"\nThis simplifies to a ratio of {non_suppressing_parts} vermilion to {suppressing_parts} wild-type.")


solve_genetics_cross()