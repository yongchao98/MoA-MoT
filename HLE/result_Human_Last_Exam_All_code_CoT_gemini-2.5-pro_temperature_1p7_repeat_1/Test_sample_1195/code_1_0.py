import math

def solve_genetics_problem():
    """
    Calculates the F2 phenotypic ratio for the given Drosophila cross.
    """
    # Step 1 & 2: Define F1 genotypes for the F1 cross
    # P Cross: Female (XvXv; su-v/su-v) x Male (XvY; Su-v/Su-v)
    # F1 Generation: All offspring are heterozygous for the suppressor gene (Su-v/su-v)
    # and have the vermilion allele (v), making them phenotypically vermilion.
    # We will cross an F1 female with an F1 male.
    f1_female = {"sex_linked": ("Xv", "Xv"), "autosomal": ("Su-v", "su-v")}
    f1_male = {"sex_linked": ("Xv", "Y"), "autosomal": ("Su-v", "su-v")}

    # Step 3: Generate all possible gametes from the F1 generation
    # F1 Female Gametes
    female_gametes = []
    for auto_allele in f1_female["autosomal"]:
        female_gametes.append(("Xv", auto_allele))

    # F1 Male Gametes
    male_gametes = []
    for sex_allele in f1_male["sex_linked"]:
        for auto_allele in f1_male["autosomal"]:
            male_gametes.append((sex_allele, auto_allele))

    # Step 4: Generate F2 generation and analyze phenotypes
    phenotype_counts = {"wild-type": 0, "vermilion": 0}

    for f_gamete in female_gametes:
        for m_gamete in male_gametes:
            # Determine autosomal genotype of the F2 offspring
            autosomal_geno = tuple(sorted((f_gamete[1], m_gamete[1])))

            # Determine phenotype based on the autosomal genotype
            # If the genotype is homozygous recessive 'su-v/su-v', the phenotype is wild-type (suppressed).
            # Otherwise, the phenotype is vermilion.
            if autosomal_geno == ("su-v", "su-v"):
                phenotype_counts["wild-type"] += 1
            else:
                phenotype_counts["vermilion"] += 1

    # Step 5: Calculate and print the final ratio
    wild_type_count = phenotype_counts["wild-type"]
    vermilion_count = phenotype_counts["vermilion"]
    total_count = wild_type_count + vermilion_count

    print("F2 Phenotypic Ratio Calculation:")
    print("---------------------------------")
    print(f"Number of vermilion offspring: {vermilion_count}")
    print(f"Number of wild-type offspring: {wild_type_count}")
    print(f"Total offspring combinations: {total_count}")
    print("---------------------------------")
    print("The final equation for the phenotypic ratio is:")
    print(f"{vermilion_count}/{total_count} vermilion : {wild_type_count}/{total_count} wild-type")
    
    # Simplify the fraction for the final answer
    common_divisor = math.gcd(vermilion_count, wild_type_count)
    simplified_vermilion = vermilion_count // common_divisor
    simplified_wild = wild_type_count // common_divisor
    
    print(f"This simplifies to the ratio: {simplified_vermilion} vermilion : {simplified_wild} wild-type")
    print(f"In fractions, this is {simplified_vermilion}/{simplified_vermilion + simplified_wild} vermilion : {simplified_wild}/{simplified_vermilion + simplified_wild} wild-type.")


solve_genetics_problem()