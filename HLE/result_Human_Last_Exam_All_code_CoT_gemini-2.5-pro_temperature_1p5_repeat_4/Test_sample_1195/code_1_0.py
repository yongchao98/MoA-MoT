def calculate_f2_ratio():
    """
    Calculates the F2 phenotypic ratio for the given Drosophila cross.
    """
    print("This genetics problem describes recessive epistasis.")
    print("The P cross is (Female) X^vX^v; su-v/su-v x (Male) X^vY; su-v+/su-v+.")
    print("All F1 offspring are X^vX^v; su-v+/su-v (females) or X^vY; su-v+/su-v (males).")
    print("Since all F1s have the vermilion allele (v) and are heterozygous for the suppressor (su-v+/su-v), they all have vermilion eyes.")
    print("\nThe F2 phenotype is determined by the cross of the F1s for the suppressor gene: su-v+/su-v x su-v+/su-v.")
    
    # F2 genotypic ratio from a monohybrid cross is 1:2:1
    homozygous_dominant = 1  # su-v+/su-v+ -> Vermilion
    heterozygous = 2         # su-v+/su-v  -> Vermilion
    homozygous_recessive = 1 # su-v/su-v   -> Wild-type (suppressed)
    total_parts = homozygous_dominant + heterozygous + homozygous_recessive

    print(f"\nThe resulting F2 genotypic ratio for the su-v gene is {homozygous_dominant}:{heterozygous}:{homozygous_recessive}.")
    
    # Calculate phenotypes
    vermilion_parts = homozygous_dominant + heterozygous
    wild_type_parts = homozygous_recessive

    print("\nPhenotype Calculation:")
    print(f"The genotypes {homozygous_dominant} (su-v+/su-v+) and {heterozygous} (su-v+/su-v) result in the vermilion phenotype.")
    print(f"The total parts for vermilion are: {homozygous_dominant} + {heterozygous} = {vermilion_parts}")

    print(f"\nThe genotype {homozygous_recessive} (su-v/su-v) results in the wild-type phenotype due to suppression.")
    
    print("\n--- Final F2 Phenotypic Ratio ---")
    print(f"The final ratio is {vermilion_parts} vermilion to {wild_type_parts} wild-type.")
    print(f"Expressed as fractions, this is: {vermilion_parts}/{total_parts} vermilion : {wild_type_parts}/{total_parts} wild-type.")

calculate_f2_ratio()