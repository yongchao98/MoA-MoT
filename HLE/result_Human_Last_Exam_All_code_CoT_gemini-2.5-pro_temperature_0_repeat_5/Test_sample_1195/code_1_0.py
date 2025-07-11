def calculate_f2_ratio():
    """
    Calculates the F2 phenotypic ratio for the given Drosophila cross.
    
    This problem involves an X-linked gene (vermilion, v) and an autosomal
    suppressor gene (su-v).

    P Cross: Female (XvXv; su-v/su-v) x Male (XvY; su-v+/su-v+)
    F1 Generation: All offspring are vermilion (XvXv; su-v+/su-v or XvY; su-v+/su-v).

    F2 Cross (F1 intercross):
    The eye color phenotype in the F2 generation depends solely on the autosomal
    suppressor gene, as all F2 flies inherit the vermilion 'v' allele.
    The cross for the suppressor gene is a standard monohybrid cross:
    su-v+/su-v x su-v+/su-v
    """

    # In a monohybrid cross (e.g., Aa x Aa), the phenotypic ratio is 3:1
    # if one allele is dominant.
    # Genotypes: 1 AA : 2 Aa : 1 aa
    # Phenotypes: 3 (Dominant) : 1 (Recessive)

    # In our case:
    # 'su-v+' is the dominant allele (no suppression).
    # 'su-v' is the recessive allele (suppression occurs).

    # Total parts in the ratio
    total_parts = 4

    # Parts with at least one dominant 'su-v+' allele.
    # These flies will have the vermilion phenotype expressed.
    vermilion_parts = 3

    # Parts with two recessive 'su-v' alleles.
    # The vermilion phenotype is suppressed, resulting in wild-type eyes.
    wild_type_parts = 1

    print("F2 Phenotypic Ratio Calculation:")
    print("The F2 phenotype is determined by the monohybrid cross of the suppressor gene: su-v+/su-v x su-v+/su-v.")
    print(f" - {vermilion_parts}/{total_parts} of the offspring will have at least one 'su-v+' allele, resulting in vermilion eyes.")
    print(f" - {wild_type_parts}/{total_parts} of the offspring will be 'su-v/su-v', resulting in suppressed (wild-type) eyes.")
    print("\nFinal F2 Phenotypic Ratio:")
    print(f"{vermilion_parts}/{total_parts} vermilion : {wild_type_parts}/{total_parts} wild-type")

calculate_f2_ratio()
<<<B>>>