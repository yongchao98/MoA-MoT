def get_reaction_product():
    """
    This function determines and prints the major product of the reaction between
    (1S,2R)-1-bromo-2-methylcyclohexane and potassium tert-butoxide.
    
    The reaction is an E2 elimination. Key points:
    1.  Substrate: (1S,2R)-1-bromo-2-methylcyclohexane is a trans-disubstituted cyclohexane.
    2.  Reagent: Potassium tert-butoxide is a strong, bulky base, favoring elimination.
    3.  Mechanism: The E2 reaction requires an anti-periplanar (diaxial) arrangement
        of the leaving group (Br) and the proton to be eliminated.
    4.  Reactive Conformer: The reaction must proceed through the less stable conformer where
        both the Br at C1 and the methyl group at C2 are in axial positions.
    5.  Regioselectivity: In this axial conformation:
        - The hydrogen at C2 is equatorial, so it cannot be eliminated. The Zaitsev
          product (1-methylcyclohexene) cannot form.
        - The hydrogen at C6 has an axial proton, which is anti-periplanar to the
          axial bromine. Elimination occurs here.
    6.  Product: Elimination of H from C6 and Br from C1 results in a double bond
        between C1 and C6.
    """
    
    # The major product of the reaction.
    # The double bond is formed between C1 and C6, with the methyl group remaining on C3
    # relative to the new double bond.
    product_name = "3-methylcyclohexene"
    
    print(f"The name of the major product is: {product_name}")

get_reaction_product()