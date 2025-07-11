def explain_regioselectivity():
    """
    Provides a detailed explanation for the observed regioselectivity in the given reactions.
    """
    explanation = """
The difference in regioselectivity between the reactions arises from the structural and electronic properties of the intermediate Münchnone (a 1,3-dipole) formed in each case.

The overall reaction is a [3+2] cycloaddition. In all cases, the N-acyl amino acid reacts with acetic anhydride (Ac2O) and triethylamine (TEA) at 65°C for 3 hours to form a Münchnone intermediate. This dipole then reacts with the dipolarophile, methyl propiolate. The regioselectivity of this final cycloaddition step depends on the symmetry of the Münchnone.

Contribution of Reactions 1 and 2:
In Reactions 1 and 2, the starting material is N-acetyl-N-methyl-alanine. The resulting Münchnone has methyl groups at both ends of its dipole, specifically at the C2 and C4 positions. These two positions are sterically and electronically very similar, creating a 'pseudo-symmetric' dipole. Consequently, there is no significant preference for how the asymmetric methyl propiolate alkyne adds across the dipole.

Reaction 2 provides the definitive evidence for this lack of selectivity. The 13C isotope labeling experiment shows that starting with a label on either the acetyl methyl or the alanine methyl results in a complete scrambling, producing a 1:1 mixture of labeled regioisomers of Product A (C9H13NO2). This demonstrates that the cycloaddition is entirely non-regioselective.

Reason for Regioselectivity in Reaction 3:
In Reaction 3, the starting material is N-acetyl-proline. The Münchnone formed from this is highly asymmetric. While position C2 still has a methyl group, position C4 is a bridgehead carbon that is part of the fused, five-membered pyrrolidine ring. This makes the C4 position much more sterically hindered and electronically distinct from C2. This structural asymmetry dictates the outcome of the cycloaddition. To minimize steric clash, the larger methoxycarbonyl group (-COOCH3) of methyl propiolate adds to the less hindered C2 position of the Münchnone. This leads to the highly regioselective formation of a single product, Product B (C10H13NO2), with no evidence of the other possible isomer.

In conclusion, Reactions 1 and 2 demonstrate that a symmetrically substituted Münchnone intermediate leads to a lack of regioselectivity. The high regioselectivity in Reaction 3 is the direct result of the significant steric and electronic asymmetry imposed on the Münchnone intermediate by the rigid, fused ring of the proline starting material.
"""
    print(explanation)

explain_regioselectivity()