import textwrap

def explain_experimental_control():
    """
    Explains the rationale for the experimental control in the given scenario.
    """

    concentration = 500  # mM GalNAc

    explanation = f"""
    1.  The Goal: The experiment aims to show that an antibody's binding to MUC1 depends on a specific sugar (Tn antigen). An inhibition assay using a high concentration of free sugar ({concentration} mM GalNAc) is performed.

    2.  The Potential Problem: A high concentration like {concentration} mM GalNAc can be toxic or cause osmotic stress to cells. This stress could cause the cells to internalize or shed their surface proteins, including the MUC1-FLAG construct being studied.

    3.  The Consequence: If the MUC1 protein is no longer on the surface, the antibody signal will decrease. This would be an experimental artifact, not the specific competitive inhibition the researcher wants to prove. The conclusion would be invalid.

    4.  The Solution (The Control): An anti-FLAG antibody binds to the protein backbone of MUC1, not the sugar. By using it, we can measure the total amount of MUC1 on the cell surface. If the anti-FLAG signal remains constant between the control (PBS) and the {concentration} mM GalNAc treated cells, it confirms that the MUC1 surface expression has not changed. This validates the experiment's main conclusion.

    5.  The Timing: The anti-FLAG antibody binds directly to the target protein on the cell, just like the main antibody of interest. Therefore, it is a primary antibody and must be added during the primary antibody incubation step to serve as a proper control for that step.

    Conclusion: Based on this reasoning, the anti-FLAG antibody is essential. It should be added with the primary antibodies to verify that the {concentration} mM GalNAc treatment has not altered the surface expression of MUC1.
    
    This corresponds to Choice C.
    """

    print(textwrap.dedent(explanation))

explain_experimental_control()