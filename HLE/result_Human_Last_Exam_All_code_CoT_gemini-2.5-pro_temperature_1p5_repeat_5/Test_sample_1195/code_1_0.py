import math

def solve_genetics_problem():
    """
    Calculates the F2 phenotypic ratio for the given Drosophila cross.
    """
    # Plan Step 5 & 6: Calculate and print the phenotypic ratio.
    # The F1 cross for the autosomal suppressor gene is su-v+/su-v x su-v+/su-v.
    # This gives a standard Mendelian genotypic ratio:
    # 1 part su-v+/su-v+
    # 2 parts su-v+/su-v
    # 1 part su-v/su-v
    
    total_parts = 1 + 2 + 1
    
    # Phenotype determination:
    # - If genotype is su-v+/su-v+ or su-v+/su-v, the vermilion color is expressed.
    # - If genotype is su-v/su-v, the suppressor restores wild-type color.
    
    vermilion_parts = 1 + 2  # Corresponds to su-v+/_ genotypes
    wild_type_parts = 1      # Corresponds to su-v/su-v genotype

    # To get the final ratio as fractions
    vermilion_fraction_numerator = vermilion_parts
    wild_type_fraction_numerator = wild_type_parts
    fraction_denominator = total_parts
    
    # Simplify the fraction if possible
    common_divisor_vermilion = math.gcd(vermilion_fraction_numerator, fraction_denominator)
    common_divisor_wild_type = math.gcd(wild_type_fraction_numerator, fraction_denominator)
    
    simplified_v_num = vermilion_fraction_numerator // common_divisor_vermilion
    simplified_v_den = fraction_denominator // common_divisor_vermilion
    
    simplified_wt_num = wild_type_fraction_numerator // common_divisor_wild_type
    simplified_wt_den = fraction_denominator // common_divisor_wild_type
    
    print("The F2 phenotypic ratio is determined by the autosomal suppressor gene (su-v).")
    print("The cross is between two heterozygotes (su-v+/su-v), resulting in a 3:1 ratio.")
    print("\n--- Final Equation ---")
    print(f"Vermilion eye color = {vermilion_parts} out of {total_parts} = {simplified_v_num}/{simplified_v_den}")
    print(f"Wild-type eye color = {wild_type_parts} out of {total_parts} = {simplified_wt_num}/{simplified_wt_den}")
    print("\nResulting Ratio:")
    print(f"{simplified_v_num}/{simplified_v_den} vermilion : {simplified_wt_num}/{simplified_wt_den} wild-type")

solve_genetics_problem()