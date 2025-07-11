def explain_regioselectivity():
    """
    This function prints the explanation for the observed regioselectivity in the chemical reactions.
    """
    explanation = """
The reason for the high regioselectivity in Reaction 3 is the steric difference between the two ends of the 1,3-dipole intermediate, a feature which is absent in the intermediates of Reactions 1 and 2.

Here is a step-by-step rationale:

1.  **Reaction Mechanism:** All three reactions proceed by converting the N-acyl amino acid starting material into a mesoionic 1,3-dipole called a münchnone. This is achieved with acetic anhydride (Ac2O) and triethylamine (TEA). This münchnone then reacts with methyl propiolate in a [3+2] cycloaddition. A subsequent loss of CO2 yields the final substituted pyrrole product.

2.  **Control Experiment (Reactions 1 & 2):** Reaction 2, using isotopically labeled starting materials, is the key experiment. It shows that when N-acetyl-N-methyl-alanine is used, the resulting münchnone has a methyl group at both its C2 and C4 positions. Because these substituents are nearly identical in size and electronic influence, the dipole is symmetric. This leads to a complete lack of regioselectivity in the cycloaddition, producing a 1:1 ratio of the two possible regioisomeric products (collectively called Product A, with formula C9H13NO2). This result demonstrates that the selectivity is governed by the substituents on the münchnone ring.

3.  **Asymmetric Dipole (Reaction 3):** In Reaction 3, the starting material is N-acetyl-proline. The münchnone formed from proline is asymmetric. While the C4 position has a methyl group, the C2 position is a bridgehead carbon that is part of the fused 5-membered pyrrolidine ring.

4.  **Steric-Controlled Selectivity:** This structural asymmetry is the reason for the high regioselectivity. The fused ring at C2 creates a much more sterically hindered environment than the methyl group at C4. During the cycloaddition, the transition state of lowest energy will be the one that minimizes steric repulsion. Therefore, the larger substituent on the alkyne (the -COOCH3 group) adds to the less-hindered C4 position of the münchnone, and the smaller hydrogen atom adds to the more-hindered C2 position. This strong steric bias results in the highly regioselective formation of a single product, Product B (with formula C10H13NO2).

In essence, Reactions 1 and 2 establish a baseline, proving that a symmetric dipole gives no selectivity (a 1:1 product ratio). Reaction 3 is selective because the proline-derived dipole is asymmetric, and steric hindrance directs the cycloaddition to a single outcome.
"""
    print(explanation)

explain_regioselectivity()