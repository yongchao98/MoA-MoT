def solve_genetics_problem():
    """
    Calculates the F2 phenotypic ratio for a Drosophila cross involving
    X-linked vermilion (v) and autosomal suppressor of vermilion (su-v).
    """

    # The problem simplifies to a monohybrid cross for the suppressor gene (su-v),
    # as all flies in the F2 generation will have the genetic basis for vermilion eyes.
    # The F1 generation cross for this gene is: su-v+/su-v x su-v+/su-v.
    
    # We analyze the F2 generation in terms of 4 total parts, based on Mendelian ratios.
    total_parts = 4
    
    # F2 genotypic ratio is 1:2:1
    # 1 part su-v+/su-v+
    # 2 parts su-v+/su-v
    # 1 part su-v/su-v
    genotype_ss_plus_homozygous = 1
    genotype_ss_plus_heterozygous = 2
    genotype_ss_homozygous = 1
    
    # Map these genotypes to phenotypes:
    # - The presence of at least one dominant su-v+ allele means no suppression occurs,
    #   and the vermilion phenotype is expressed.
    # - The homozygous recessive su-v/su-v genotype suppresses vermilion,
    #   resulting in a wild-type phenotype.
    
    vermilion_count = genotype_ss_plus_homozygous + genotype_ss_plus_heterozygous
    wild_type_count = genotype_ss_homozygous

    print("Step 1: The F1 cross for the suppressor gene is su-v+/su-v x su-v+/su-v.")
    print(f"Step 2: The F2 genotypic ratio is {genotype_ss_plus_homozygous} (su-v+/su-v+) : {genotype_ss_plus_heterozygous} (su-v+/su-v) : {genotype_ss_homozygous} (su-v/su-v).")
    print("\nStep 3: Calculating phenotypes from genotypes:")
    print(f" - Vermilion (no suppression) = {genotype_ss_plus_homozygous} part(s) + {genotype_ss_plus_heterozygous} part(s) = {vermilion_count} parts.")
    print(f" - Wild-type (suppression) = {genotype_ss_homozygous} part(s).")
    
    print("\nFinal Phenotypic Ratio Equation:")
    # The final equation expresses the ratio as fractions of the total.
    print(f"{vermilion_count}/{total_parts} vermilion : {wild_type_count}/{total_parts} wild-type")

solve_genetics_problem()