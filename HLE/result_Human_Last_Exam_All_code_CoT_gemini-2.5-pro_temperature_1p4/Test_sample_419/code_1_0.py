import textwrap

def explain_experimental_control():
    """
    Explains the reasoning for the correct experimental control in the given scenario.
    """

    reasoning = """
    1. The goal of the experiment is to show that free GalNAc sugar competitively inhibits the antibody's binding to a glycan on the MUC1 protein.

    2. A potential confounding variable is that a high concentration (500 mM) of GalNAc might have other, "off-target" effects on the cells. Specifically, it could be toxic or cause cellular stress, leading to the internalization or shedding of the MUC1 protein from the cell surface.

    3. If the MUC1 protein is no longer on the surface, the antibody cannot bind. This would lead to a lower signal, but the reason would be the absence of the target, not the inhibition of the binding event. This would invalidate the conclusion.

    4. To control for this, it is essential to measure the amount of MUC1 protein on the surface, independent of the sugar-dependent antibody. The MUC1 construct has a FLAG tag for this purpose.

    5. An anti-FLAG antibody binds directly to the FLAG tag on the MUC1 protein. This binding is not dependent on the sugar. By using the anti-FLAG antibody, one can confirm that the total amount of MUC1 on the cell surface is the same in both the PBS (control) and GalNAc (treatment) conditions.

    6. The anti-FLAG antibody is a primary antibody, just like the main antibody being tested. In a standard immunoassay, all primary antibodies are added at the same time (the primary incubation step).

    7. Therefore, the anti-FLAG antibody should be added with the primary antibodies to verify that the GalNAc treatment has not altered the surface expression of the MUC1 protein. This logic directly corresponds to answer choice C.
    """
    print(textwrap.dedent(reasoning).strip())
    print("\nBased on this reasoning, the correct answer is C.")
    print("<<<C>>>")

explain_experimental_control()