def explain_regioselectivity():
    """
    Provides a detailed explanation for the observed regioselectivity in the chemical reactions.
    """
    explanation = """
The reason for the high regioselectivity in Reaction 3 is the inherent asymmetry of its münchnone intermediate, with the results from Reactions 1 and 2 providing the critical evidence to support this conclusion.

1.  **Reaction 3 (High Regioselectivity):** The starting material, N-acetylproline, is cyclic. Upon treatment with acetic anhydride and TEA, it forms a rigid, bicyclic münchnone intermediate. In this intermediate, the two reacting carbons of the 1,3-dipole, C2 and C5, are structurally non-equivalent. C2 bears a methyl group, while C5 is a tertiary bridgehead carbon that is part of the original proline ring. This permanent steric and electronic asymmetry forces the unsymmetrical alkyne (methyl propiolate) to add in a specific orientation, leading to the formation of a single regioisomer, Product B.

2.  **Reactions 1 & 2 (The Control Experiment):** The starting material for these reactions is the acyclic N-acetyl-N-methyl-alanine. This forms a münchnone intermediate where both C2 and C5 have methyl substituents. While these methyl groups originate from different parts of the starting material (acetyl vs. alanine), the crucial insight comes from Reaction 2.

3.  **The Role of the Isotope Labeling (Reaction 2):** The 13C labeling experiment demonstrates that the labels placed on the acetyl-methyl and the alanine-methyl are completely scrambled in the final product, appearing in a 1:1 ratio at both methyl positions of the pyrrole ring. This definitively proves that the münchnone intermediate is effectively **symmetrical**; the C2 and C5 positions become chemically indistinguishable during the reaction.

**Conclusion:**

The regioselectivity in Reaction 3 is caused by the **asymmetry of the proline-derived münchnone**. Reactions 1 and 2 establish a baseline for a system with a **symmetrical münchnone**. The isotope scrambling in Reaction 2 proves this symmetry. Because the dipole in Reaction 1 is symmetrical, it lacks the structural bias needed to control the orientation of the incoming unsymmetrical alkyne. This contrasts sharply with Reaction 3, where the dipole's built-in asymmetry acts as a "steering" group, ensuring the reaction proceeds with high regioselectivity.
"""
    print(explanation.strip())

explain_regioselectivity()