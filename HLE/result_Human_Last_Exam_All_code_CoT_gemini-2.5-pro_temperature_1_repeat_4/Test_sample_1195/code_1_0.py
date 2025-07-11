def solve_drosophila_cross():
    """
    Calculates the F2 phenotypic ratio for a cross involving an X-linked
    gene (vermilion) and an autosomal suppressor gene (su-v).
    The code follows the logical steps derived from Mendelian genetics.
    """

    # Define parental genotypes
    p_female_genotype = "X(v)X(v); su-v/su-v"
    p_male_genotype = "X(v)Y; su-v+/su-v+"

    # F1 genotypes are determined by the parental cross
    f1_female_genotype = "X(v)X(v); su-v+/su-v"
    f1_male_genotype = "X(v)Y; su-v+/su-v"

    print("--- Problem Analysis ---")
    print(f"Parental Female (P): {p_female_genotype}")
    print(f"Parental Male (P):   {p_male_genotype}")
    print(f"F1 Generation Genotype (Male and Female): ...; su-v+/su-v")
    print("\nSince all F1 flies have a vermilion genotype (Xv) and are heterozygous for the suppressor,")
    print("the F2 phenotype depends only on the segregation of the su-v gene.")
    print("The F2 cross for the suppressor gene is: su-v+/su-v x su-v+/su-v\n")

    # For a standard monohybrid cross (Ss x Ss), the offspring genotypic ratio is 1:2:1.
    # We represent this ratio using integers for clarity.
    # 1 part su-v+/su-v+
    # 2 parts su-v+/su-v
    # 1 part su-v/su-v
    prop_dominant_homozygous = 1
    prop_heterozygous = 2
    prop_recessive_homozygous = 1
    total_parts = prop_dominant_homozygous + prop_heterozygous + prop_recessive_homozygous

    # Determine phenotype based on genotype
    # su-v+/su-v+ and su-v+/su-v do not suppress -> vermilion
    # su-v/su-v suppresses -> wild-type
    vermilion_parts = prop_dominant_homozygous + prop_heterozygous
    wild_type_parts = prop_recessive_homozygous

    print("--- F2 Phenotypic Ratio Calculation ---")
    print(f"The genotypic ratio for the su-v gene is {prop_dominant_homozygous} (su-v+/su-v+) : {prop_heterozygous} (su-v+/su-v) : {prop_recessive_homozygous} (su-v/su-v).")

    print(f"\n- Phenotype: vermilion (no suppression)")
    print(f"  - Genotypes: su-v+/su-v+ and su-v+/su-v")
    print(f"  - Contribution to ratio: {prop_dominant_homozygous} + {prop_heterozygous} = {vermilion_parts}")

    print(f"\n- Phenotype: wild-type (suppression)")
    print(f"  - Genotype: su-v/su-v")
    print(f"  - Contribution to ratio: {prop_recessive_homozygous}")

    print("\n--- Final Equation and Ratio ---")
    print(f"Out of {total_parts} total theoretical offspring in the F2 generation:")
    print(f"The equation for the phenotypes is: {vermilion_parts} vermilion + {wild_type_parts} wild-type = {total_parts} total")
    print(f"This gives a phenotypic ratio of {vermilion_parts}/{total_parts} vermilion to {wild_type_parts}/{total_parts} wild-type.")
    print(f"The final simplified ratio is 3/4 vermilion : 1/4 wild-type.")


if __name__ == '__main__':
    solve_drosophila_cross()