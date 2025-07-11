import math

def calculate_f2_phenotypic_ratio():
    """
    Calculates the F2 phenotypic ratio for the given Drosophila cross.

    The problem involves an X-linked gene (vermilion, v) and an autosomal
    suppressor gene (su-v). Since all flies in the relevant generations have
    the 'v' allele, the phenotype is determined by the su-v locus.

    The F1 generation is heterozygous for the suppressor gene (su-v+/su-v).
    The F2 generation results from an F1 x F1 cross (su-v+/su-v x su-v+/su-v).
    """

    # F2 genotypic ratio from a su-v+/su-v x su-v+/su-v cross is 1:2:1
    # Proportion of su-v+/su-v+
    prop_homozygous_dominant = 1
    # Proportion of su-v+/su-v
    prop_heterozygous = 2
    # Proportion of su-v/su-v (recessive suppressor is active)
    prop_homozygous_recessive = 1

    # Total parts in the ratio
    total_parts = prop_homozygous_dominant + prop_heterozygous + prop_homozygous_recessive

    # Phenotype calculation
    # Vermilion phenotype occurs when the suppressor is not active (su-v+/su-v+ or su-v+/su-v)
    vermilion_count = prop_homozygous_dominant + prop_heterozygous
    
    # Wild-type phenotype occurs when the suppressor is active (su-v/su-v)
    wild_type_count = prop_homozygous_recessive
    
    # Print the explanation and result
    print("In the F2 generation, the phenotypic ratio is determined by the autosomal suppressor gene (su-v).")
    print("The genotypic ratio for this gene is 1 (su-v+/su-v+) : 2 (su-v+/su-v) : 1 (su-v/su-v).")
    print("-" * 30)
    print(f"Number of 'parts' resulting in Vermilion phenotype (su-v+/_): {prop_homozygous_dominant} + {prop_heterozygous} = {vermilion_count}")
    print(f"Number of 'parts' resulting in Wild-type phenotype (su-v/su-v): {wild_type_count}")
    print("-" * 30)
    print("The final F2 phenotypic ratio is derived from these counts.")
    print(f"\nFinal Answer: {vermilion_count}/{total_parts} vermilion : {wild_type_count}/{total_parts} wild-type")


# Execute the function to get the answer
calculate_f2_phenotypic_ratio()