import math

def get_phenotype(sex_linked, autosomal):
    """
    Determines the eye color phenotype based on the genotype.
    - vermilion (v) is an X-linked recessive gene.
    - suppressor of vermilion (su-v) is an autosomal recessive gene.
    - The su-v/su-v genotype restores wild-type eye color in flies
      that are genetically vermilion.
    """
    # Is the fly genetically vermilion? (XvXv or XvY)
    is_vermilion_genotype = "Xv" in sex_linked
    
    # Does the fly have the homozygous recessive suppressor genotype?
    has_suppressor = autosomal == "su-v/su-v"
    
    if is_vermilion_genotype and has_suppressor:
        return "wild-type"  # Suppression restores wild-type color
    elif is_vermilion_genotype:
        return "vermilion" # No suppression, vermilion is expressed
    else:
        # This case is not reached in this specific problem
        return "wild-type"

def solve_genetics_cross():
    """
    Solves the Drosophila eye color cross problem step-by-step.
    """
    print("--- Solving the Drosophila Gene Interaction Problem ---")
    
    # Step 1 & 2: Define P generation and determine F1 generation
    print("\nStep 1: Parental (P) cross and resulting F1 generation.")
    p_female_geno = "XvXv; su-v/su-v"
    p_male_geno = "XvY; +/+"
    f1_female_geno = "XvXv; +/su-v"
    f1_male_geno = "XvY; +/su-v"
    
    print(f"  - P Cross: Female ({p_female_geno}) x Male ({p_male_geno})")
    print(f"  - All F1 offspring are heterozygous for the suppressor gene.")
    print(f"  - F1 Female genotype: {f1_female_geno}")
    print(f"  - F1 Male genotype:   {f1_male_geno}")
    print("  - F1 Phenotype: Since the suppressor allele 'su-v' is recessive, all F1 flies are vermilion.")

    # Step 3 & 4: F1 cross and gamete formation
    print("\nStep 2: F1 cross (XvXv; +/su-v  x  XvY; +/su-v) to produce F2 generation.")
    
    # Gametes from F1 female (XvXv; +/su-v) -> 1/2 Xv;+ and 1/2 Xv;su-v
    f1_female_gametes = [("Xv", "+"), ("Xv", "su-v")]
    
    # Gametes from F1 male (XvY; +/su-v) -> 1/4 Xv;+, 1/4 Xv;su-v, 1/4 Y;+, 1/4 Y;su-v
    f1_male_gametes = [("Xv", "+"), ("Xv", "su-v"), ("Y", "+"), ("Y", "su-v")]
    
    # Step 5: Calculate F2 genotypes and phenotypes
    phenotype_counts = {"wild-type": 0, "vermilion": 0}
    total_parts = len(f1_female_gametes) * len(f1_male_gametes)

    for f_gamete in f1_female_gametes:
        for m_gamete in f1_male_gametes:
            sex_linked_offspring = f_gamete[0] + m_gamete[0]
            
            # Combine autosomal alleles and sort them for consistency
            autosomal_alleles = sorted([f_gamete[1], m_gamete[1]])
            autosomal_offspring = f"{autosomal_alleles[0]}/{autosomal_alleles[1]}"
            
            phenotype = get_phenotype(sex_linked_offspring, autosomal_offspring)
            phenotype_counts[phenotype] += 1
            
    # Step 6: Calculate and display the final ratio
    print("\nStep 3: Calculating the F2 phenotypic ratio.")
    wild_type_count = phenotype_counts["wild-type"]
    vermilion_count = phenotype_counts["vermilion"]
    
    print(f"  - Out of {total_parts} possible outcomes:")
    print(f"    - Wild-type count: {wild_type_count}")
    print(f"    - Vermilion count: {vermilion_count}")
    
    common_divisor = math.gcd(wild_type_count, vermilion_count)
    simple_wild = wild_type_count // common_divisor
    simple_vermilion = vermilion_count // common_divisor
    
    print("\n--- Final Answer ---")
    print("The final phenotypic ratio in the F2 generation is:")
    print(f"{simple_vermilion} vermilion : {simple_wild} wild-type")
    
    total_simple = simple_wild + simple_vermilion
    print(f"\nThis is equivalent to the fractional ratio:")
    print(f"{simple_vermilion}/{total_simple} vermilion : {simple_wild}/{total_simple} wild-type")

# Execute the solution
solve_genetics_cross()
<<<B>>>