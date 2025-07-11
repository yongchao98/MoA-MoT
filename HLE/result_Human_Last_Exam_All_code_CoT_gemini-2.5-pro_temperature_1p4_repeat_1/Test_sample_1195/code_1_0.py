import collections

def solve_genetics_ratio():
    """
    Calculates the F2 phenotypic ratio for the given Drosophila cross.
    The key insight is that the F2 phenotype is determined solely by the
    autosomal 'su-v' gene cross, as all F2 flies are genetically vermilion.
    """
    
    # Define the alleles of the F1 parents for the suppressor gene locus.
    # The F1 cross is su-v+/su-v x su-v+/su-v.
    f1_parent_alleles = ['su-v+', 'su-v']

    # Generate all possible F2 genotypes for the su-v locus.
    f2_genotypes = []
    for allele1 in f1_parent_alleles:
        for allele2 in f1_parent_alleles:
            # Sorting alleles helps group heterozygotes like 'su-v+/su-v' and 'su-v/su-v+' together.
            genotype = '/'.join(sorted([allele1, allele2]))
            f2_genotypes.append(genotype)

    # Count the occurrences of each F2 genotype.
    genotype_counts = collections.Counter(f2_genotypes)
    
    # Determine phenotypes based on genotypes.
    # 'su-v/su-v' suppresses the vermilion trait, resulting in wild-type eyes.
    # All other genotypes result in vermilion eyes.
    phenotype_counts = collections.Counter()
    phenotype_counts['wild-type'] = genotype_counts['su-v/su-v']
    phenotype_counts['vermilion'] = genotype_counts['su-v+/su-v+'] + genotype_counts['su-v+/su-v']

    total_offspring = len(f2_genotypes)
    wild_type_num = phenotype_counts['wild-type']
    vermilion_num = phenotype_counts['vermilion']
    
    print("F2 Phenotypic Ratio Calculation:")
    print("-" * 30)
    print("The F1 cross for the autosomal suppressor gene is: su-v+/su-v  x  su-v+/su-v")
    print("\nThe resulting F2 genotypic ratio for this gene is:")
    print(f"1 su-v+/su-v+ : 2 su-v+/su-v : 1 su-v/su-v")
    print("\nPhenotype mapping:")
    print(f"  - Genotypes su-v+/su-v+ and su-v+/su-v do not suppress -> vermilion eyes.")
    print(f"  - Genotype su-v/su-v suppresses -> wild-type eyes.")
    print("\nFinal Ratio Equation:")
    print(f"The fraction of vermilion offspring is {vermilion_num}/{total_offspring}.")
    print(f"The fraction of wild-type offspring is {wild_type_num}/{total_offspring}.")
    print(f"\nTherefore, the final phenotypic ratio is {vermilion_num} vermilion : {wild_type_num} wild-type.")

solve_genetics_ratio()