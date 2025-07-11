import math
from collections import Counter

def get_phenotype(genotype):
    """Determines the phenotype based on the full genotype."""
    # The X-linked gene is always vermilion in this F2 generation (XvXv or XvY).
    # The phenotype is therefore determined by the autosomal suppressor gene.
    autosomal_genotype = genotype[1]
    
    # The suppressor is only active when homozygous recessive (su-v/su-v).
    # We sort the alleles to handle both ('su-v', 'su-v+') and ('su-v+', 'su-v') as heterozygous.
    if tuple(sorted(autosomal_genotype)) == ('su-v', 'su-v'):
        return "Wild-type"  # Suppressor is active, restoring wild-type color
    else:
        return "Vermilion"  # Suppressor is not active

def solve_genetics_cross():
    """
    Solves the Drosophila cross problem to find the F2 phenotypic ratio.
    """
    # Step 1 & 2: Define F1 genotypes based on the P cross
    # P Cross: XvXv su-v/su-v (female) x XvY su-v+/su-v+ (male)
    # F1 Genotypes are all XvXv su-v+/su-v (females) and XvY su-v+/su-v (males).
    f1_female = (("Xv", "Xv"), ("su-v+", "su-v"))
    f1_male = (("Xv", "Y"), ("su-v+", "su-v"))

    # Step 3: Determine gametes from the F1 generation
    # F1 Female (XvXv su-v+/su-v) produces two gamete types:
    female_gametes = [("Xv", "su-v+"), ("Xv", "su-v")]
    
    # F1 Male (XvY su-v+/su-v) produces four gamete types:
    male_gametes = [("Xv", "su-v+"), ("Xv", "su-v"), ("Y", "su-v+"), ("Y", "su-v")]

    # Step 4: Generate all possible F2 offspring genotypes
    f2_genotypes = []
    for f_gamete in female_gametes:
        for m_gamete in male_gametes:
            sex_alleles = (f_gamete[0], m_gamete[0])
            autosomal_alleles = (f_gamete[1], m_gamete[1])
            f2_genotypes.append((sex_alleles, autosomal_alleles))
            
    # Step 5: Determine and count F2 phenotypes
    f2_phenotypes = [get_phenotype(geno) for geno in f2_genotypes]
    phenotype_counts = Counter(f2_phenotypes)
    
    total_parts = len(f2_phenotypes)
    
    # Print the results and the equation
    print("F2 Phenotypic Ratio Calculation:\n")
    
    vermilion_count = phenotype_counts.get("Vermilion", 0)
    wild_type_count = phenotype_counts.get("Wild-type", 0)
    
    # Calculate simplified fraction for vermilion
    gcd_v = math.gcd(vermilion_count, total_parts)
    print(f"Vermilion flies = {vermilion_count} / {total_parts}")
    print(f"Simplified: {vermilion_count // gcd_v} / {total_parts // gcd_v}\n")
    
    # Calculate simplified fraction for wild-type
    gcd_w = math.gcd(wild_type_count, total_parts)
    print(f"Wild-type flies = {wild_type_count} / {total_parts}")
    print(f"Simplified: {wild_type_count // gcd_w} / {total_parts // gcd_w}\n")
    
    # Find simplified ratio (e.g., 3:1)
    common_divisor = math.gcd(vermilion_count, wild_type_count)
    
    print(f"The final phenotypic ratio is {vermilion_count // common_divisor} vermilion : {wild_type_count // common_divisor} wild-type.")


solve_genetics_cross()