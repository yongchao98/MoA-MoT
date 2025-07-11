import textwrap

def explain_experimental_control():
    """
    Explains the reasoning behind the choice of experimental control.
    """

    explanation = """
    The goal is to confirm that the observed inhibition of antibody binding is due to competition by free GalNAc, and not due to an unintended side-effect of the high GalNAc concentration.

    1.  **Potential Confounder:** A high concentration of a substance like 500 mM GalNAc could be toxic or cause stress to the cells, potentially leading them to internalize or shed surface proteins like the MUC1 construct. If this were to happen, the antibody would fail to bind because its target protein is no longer on the surface, not because of specific competition.

    2.  **The Control:** The MUC1 protein is engineered with a FLAG tag. An anti-FLAG antibody binds to this tag on the protein backbone, completely independent of the sugar (TN antigen). Therefore, the anti-FLAG antibody serves as a perfect control to measure the amount of MUC1 protein present on the cell surface.

    3.  **Timing:** The anti-FLAG antibody is a 'primary' antibody because it binds directly to its target on the cell. It must be added during the primary antibody incubation step to accurately assess the MUC1 levels under the experimental conditions.

    4.  **Conclusion:** By comparing the anti-FLAG binding in the presence of PBS versus 500 mM GalNAc, you can verify that the surface level of the MUC1 protein remains constant. This ensures that any decrease in the primary antibody's binding is a direct result of competition from the free sugar. This makes C the correct choice.
    """

    print(textwrap.dedent(explanation).strip())
    print("\nCorrect Answer Choice:")
    print("C. Anti-flag should be added with the primary antibodies. Its included to verify GalNAc has not altered surface expression of MUC1")

explain_experimental_control()