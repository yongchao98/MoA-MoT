def solve_genetics_problem():
    """
    Calculates and explains the F2 phenotypic ratio for the given Drosophila cross.
    """
    print("This script calculates the F2 phenotypic ratio for a dihybrid cross in Drosophila.")
    print("The cross involves an X-linked gene (vermilion, v) and an autosomal suppressor gene (su-v).\n")

    print("Step 1: Parental (P) Generation Cross")
    print("  - Female: X(v)X(v); su-v/su-v (Phenotype: wild-type, due to suppression)")
    print("  - Male:   X(v)Y; su-v+/su-v+ (Phenotype: vermilion)\n")

    print("Step 2: F1 Generation")
    print("  - All F1 Females are X(v)X(v); su-v+/su-v")
    print("  - All F1 Males are   X(v)Y; su-v+/su-v")
    print("  - Since all F1 individuals have the su-v+ allele, there is no suppression.")
    print("  - Therefore, all F1 flies have vermilion eyes.\n")

    print("Step 3: F2 Generation (from F1 x F1 cross)")
    print("  - We analyze the offspring from crossing F1 females and F1 males.")
    print("  - The possible F2 genotypes are calculated below:\n")

    # F2 genotypes and their proportions (out of 8 total parts)
    f2_genotypes = {
        "X(v)X(v); su-v+/su-v+": 1,
        "X(v)X(v); su-v+/su-v": 2,
        "X(v)X(v); su-v/su-v": 1,
        "X(v)Y; su-v+/su-v+": 1,
        "X(v)Y; su-v+/su-v": 2,
        "X(v)Y; su-v/su-v": 1,
    }

    wild_type_count = 0
    vermilion_count = 0

    print("Step 4: Determine F2 Phenotypes and Final Ratio\n")
    print("Phenotype determination for each F2 genotype:")
    for genotype, count in f2_genotypes.items():
        # The su-v/su-v genotype suppresses the vermilion phenotype, resulting in wild-type eyes.
        if "su-v/su-v" in genotype:
            phenotype = "Wild-type (suppressed)"
            wild_type_count += count
        else:
            phenotype = "Vermilion"
            vermilion_count += count
        print(f"  - Genotype: {genotype} (Count: {count}/8) -> Phenotype: {phenotype}")

    total_parts = wild_type_count + vermilion_count

    print("\n--- Final Calculation ---")
    print(f"Total parts resulting in Wild-type eyes = {wild_type_count}")
    print(f"Total parts resulting in Vermilion eyes = {vermilion_count}")
    print("\nThe final equation for the phenotypic ratio (Vermilion : Wild-type) is:")
    print(f"{vermilion_count} : {wild_type_count}")
    
    # Simplify the ratio
    common_divisor = wild_type_count 
    simple_vermilion = vermilion_count // common_divisor
    simple_wild_type = wild_type_count // common_divisor
    
    print(f"This simplifies to the ratio {simple_vermilion} : {simple_wild_type}.")
    
    print(f"\nTherefore, the expected phenotypic ratio in the F2 generation is:")
    print(f"{vermilion_count}/{total_parts} vermilion : {wild_type_count}/{total_parts} wild-type")
    print(f"Which simplifies to 3/4 vermilion : 1/4 wild-type.")

solve_genetics_problem()