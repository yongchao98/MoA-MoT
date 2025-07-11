import math

def solve_genetics_problem():
    """
    This function solves the Drosophila genetics problem by simulating the cross
    and calculating the F2 phenotypic ratio.
    """
    
    # Step 1: Define F1 genotypes based on the P-generation cross.
    # P-Cross: X(v)X(v); su-v/su-v (female) x X(v)Y; su-v+/su-v+ (male)
    # This results in all F1 females being X(v)X(v); su-v+/su-v
    # and all F1 males being X(v)Y; su-v+/su-v
    
    # Step 2: Simulate the F1 cross: X(v)X(v); su-v+/su-v x X(v)Y; su-v+/su-v
    # We can analyze the autosomal gene cross separately.
    # Cross: su-v+/su-v x su-v+/su-v
    # Possible offspring genotypes for the suppressor gene:
    # 1/4 su-v+/su-v+
    # 2/4 su-v+/su-v
    # 1/4 su-v/su-v
    
    # Step 3: Count phenotypes.
    # All F2 offspring will have the X(v) allele.
    # The phenotype is determined by the suppressor gene.
    # - If genotype is su-v/su-v, vermilion is suppressed -> wild-type eyes.
    # - If genotype is su-v+/_ (su-v+/su-v+ or su-v+/su-v), vermilion is expressed.
    
    # Fractions of genotypes
    prob_suppressor_homozygous = 1
    prob_not_suppressor = 3 # (1 + 2)
    total_parts = 4 # (1 + 2 + 1)

    # Phenotype counts based on the 4 possible outcomes of the autosomal cross
    wild_type_count = prob_suppressor_homozygous
    vermilion_count = prob_not_suppressor
    
    # Simplify the ratio
    common_divisor = math.gcd(wild_type_count, vermilion_count)
    
    print("Step-by-step derivation of the F2 phenotypic ratio:")
    print("-" * 60)
    print("1. Parental (P) & F1 Generation Analysis:")
    print("  - P-Cross: X(v)X(v); su-v/su-v   x   X(v)Y; su-v+/su-v+")
    print("  - All F1 offspring are heterozygous for the suppressor gene (su-v+/su-v).")
    print("-" * 60)
    print("2. F1 Cross to Produce F2 Generation:")
    print("  - F1-Cross: X(v)X(v); su-v+/su-v   x   X(v)Y; su-v+/su-v")
    print("  - Since all F1 parents have the X(v) allele, all F2 offspring will have the genetic potential for vermilion eyes.")
    print("  - The final phenotype depends entirely on the autosomal suppressor gene.")
    print("-" * 60)
    print("3. Autosomal Cross Analysis (su-v+/su-v  x  su-v+/su-v):")
    print(f"  - The probability of a non-suppressing genotype (su-v+/_): {prob_not_suppressor}/{total_parts}")
    print(f"  - The probability of a suppressing genotype (su-v/su-v): {prob_suppressor_homozygous}/{total_parts}")
    print("-" * 60)
    print("4. Final F2 Phenotypic Ratio Calculation:")
    print("  - Offspring with a non-suppressing genotype will have VERMILION eyes.")
    print("  - Offspring with the su-v/su-v genotype will have their vermilion suppressed, resulting in WILD-TYPE eyes.")
    print("\nFinal Equation:")
    print(f"The expected ratio is {vermilion_count} vermilion : {wild_type_count} wild-type.")
    print(f"Simplified, the ratio is {vermilion_count // common_divisor}/4 vermilion : {wild_type_count // common_divisor}/4 wild-type.")

solve_genetics_problem()
<<<B>>>