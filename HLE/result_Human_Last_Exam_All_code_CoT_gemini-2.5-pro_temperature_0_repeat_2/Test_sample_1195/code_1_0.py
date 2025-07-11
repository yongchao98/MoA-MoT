from fractions import Fraction

def solve_genetics_cross():
    """
    Calculates the F2 phenotypic ratio for the given Drosophila cross.
    
    Gene Notation:
    - X-linked eye color: 'v' (vermilion, recessive)
    - Autosomal suppressor: 'su_v' (suppressor, recessive), 'su_v_plus' (wild-type, dominant)
    """

    # Step 1: Define F1 genotypes based on the problem description
    # P Cross: (XvXv; su_v/su_v) x (XvY; su_v_plus/su_v_plus)
    # F1 Female: XvXv; su_v_plus/su_v
    # F1 Male: XvY; su_v_plus/su_v
    
    # Step 2: Determine the gametes produced by the F1 generation
    # F1 Female (XvXv; su_v_plus/su_v) produces:
    f1_female_gametes = {
        ('Xv', 'su_v_plus'): Fraction(1, 2),
        ('Xv', 'su_v'): Fraction(1, 2)
    }
    
    # F1 Male (XvY; su_v_plus/su_v) produces:
    f1_male_gametes = {
        ('Xv', 'su_v_plus'): Fraction(1, 4),
        ('Xv', 'su_v'): Fraction(1, 4),
        ('Y', 'su_v_plus'): Fraction(1, 4),
        ('Y', 'su_v'): Fraction(1, 4)
    }

    # Step 3: Calculate F2 generation phenotypes and their probabilities
    phenotype_probabilities = {'Wild-type': Fraction(0), 'Vermilion': Fraction(0)}

    # Iterate through all possible gamete combinations to form F2 offspring
    for f_gamete, f_prob in f1_female_gametes.items():
        for m_gamete, m_prob in f1_male_gametes.items():
            
            f_auto_allele = f_gamete[1]
            m_auto_allele = m_gamete[1]
            
            # Determine the phenotype based on the suppressor gene
            # If the genotype is su_v/su_v, the phenotype is wild-type (suppressed).
            # Otherwise, the vermilion trait is expressed.
            if f_auto_allele == 'su_v' and m_auto_allele == 'su_v':
                phenotype = 'Wild-type'
            else:
                # All F1 gametes carry the 'v' allele on the X chromosome,
                # so any non-suppressed offspring will be vermilion.
                phenotype = 'Vermilion'
            
            # Probability of this specific offspring genotype
            offspring_prob = f_prob * m_prob
            
            # Add the probability to the corresponding phenotype
            phenotype_probabilities[phenotype] += offspring_prob

    # Step 4: Print the final results
    wild_type_frac = phenotype_probabilities['Wild-type']
    vermilion_frac = phenotype_probabilities['Vermilion']

    print("F2 Phenotypic Ratio Calculation:")
    print(f"The probability of wild-type eyes (due to su-v/su-v) is {wild_type_frac.numerator}/{wild_type_frac.denominator}.")
    print(f"The probability of vermilion eyes is {vermilion_frac.numerator}/{vermilion_frac.denominator}.")
    print("\nFinal expected phenotypic ratio:")
    # The problem asks to output each number in the final equation
    print(f"{vermilion_frac.numerator}/{vermilion_frac.denominator} vermilion : {wild_type_frac.numerator}/{wild_type_frac.denominator} wild-type")

solve_genetics_cross()
<<<B>>>