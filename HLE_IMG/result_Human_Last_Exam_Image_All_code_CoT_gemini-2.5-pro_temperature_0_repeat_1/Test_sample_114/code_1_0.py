def explain_regioselectivity():
    """
    Explains the difference in regioselectivity between the given reactions.
    """
    explanation = """
The difference in regioselectivity between Reaction 3 and Reactions 1 & 2 is due to the symmetry of the key reactive intermediate, a mesoionic compound called a münchnone.

**General Mechanism:**
All three reactions proceed via a [3+2] cycloaddition. The N-acyl amino acid is first treated with acetic anhydride (Ac2O) and triethylamine (TEA) to form a münchnone. This münchnone is a 1,3-dipole that reacts with the alkyne (methyl propiolate). The resulting bicyclic adduct is unstable and rapidly loses a molecule of carbon dioxide (CO2) to form the final pyrrole product.

**Analysis of Reactions 1 and 2:**
1.  **Intermediate Formation:** The starting material, N-acetyl-N-methyl-alanine, forms a münchnone that has methyl groups at both the C2 and C4 positions of the oxazolium ring.
2.  **Symmetry:** The two ends of the 1,3-dipole (C2 and C4) are sterically and electronically very similar because both are substituted with a methyl group. This makes the reactive dipole effectively symmetric.
3.  **Lack of Selectivity:** When this symmetric münchnone reacts with the asymmetric alkyne (methyl propiolate, H-C≡C-CO2Me), there is no significant preference for which end of the dipole attacks which end of the alkyne. This lack of preference is experimentally confirmed in Reaction 2, where the isotopic labeling shows that the two possible regioisomers are formed in a 1:1 ratio.

**Analysis of Reaction 3:**
1.  **Intermediate Formation:** The starting material, N-acetyl-proline, is a cyclic amino acid. It forms a bicyclic münchnone.
2.  **Asymmetry:** In this münchnone, the C2 position is substituted with a methyl group (from the acetyl group), but the C4 position is substituted with a hydrogen atom (the original alpha-proton of proline).
3.  **High Regioselectivity:** This makes the 1,3-dipole highly asymmetric. The two ends (C2-Me vs. C4-H) are now sterically and electronically distinct. The cycloaddition is therefore controlled. The bulkier ester group (-CO2Me) of the alkyne preferentially adds to the less sterically hindered C4-H side of the münchnone. This leads to the exclusive formation of a single regioisomer, Product B.

**Conclusion:**
The results from Reactions 1 and 2 are critical because they establish a baseline: a symmetric münchnone leads to no regioselectivity (a 1:1 product ratio). This demonstrates that the high regioselectivity observed in Reaction 3 is a direct consequence of the inherent asymmetry of the proline-derived münchnone intermediate. The steric and electronic differences between the methyl-substituted C2 and the hydrogen-substituted C4 of the dipole in Reaction 3 dictate the orientation of the cycloaddition, resulting in the formation of a single product.
"""
    print(explanation)

explain_regioselectivity()