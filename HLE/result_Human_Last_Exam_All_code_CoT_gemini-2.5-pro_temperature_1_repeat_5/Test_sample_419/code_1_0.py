import textwrap

def explain_experimental_control():
    """
    Explains the reasoning behind the correct answer for the biological experiment design question.
    """
    explanation = """
    1.  **Purpose of the Anti-FLAG Antibody:** The MUC1 protein in this experiment has a FLAG tag, a small peptide sequence engineered into the protein. An anti-FLAG antibody binds directly to this tag. This binding is independent of the glycosylation state (the TN antigen) of MUC1. Therefore, the anti-FLAG antibody serves as a control to measure the total amount of MUC1 protein present on the cell surface.

    2.  **Controlling for Off-Target Effects:** The experiment uses a very high concentration of GalNAc (500 mM) as an inhibitor. Such a high concentration could have unintended, or "off-target," effects on the cells. A major concern is that it could be toxic or cause osmotic stress, leading the cells to internalize surface proteins, including the MUC1 construct.

    3.  **Verifying the Mechanism of Inhibition:** If the GalNAc treatment caused MUC1 to be removed from the cell surface, the main antibody (which recognizes the sugar on MUC1) would fail to bind. This would look like inhibition, but it would be due to the target's absence, not competitive binding. By measuring the anti-FLAG signal in both the control (PBS) and GalNAc-treated samples, we can check if the total MUC1 on the surface remains constant. If the anti-FLAG signal is unchanged, we can be confident that the MUC1 is still there and that the decreased signal from the main antibody is truly due to competitive inhibition at the binding site.

    4.  **Correct Timing of Antibody Addition:** The anti-FLAG antibody is a primary antibody because it binds directly to its target on the protein. In a typical immunoassay, all primary antibodies are added at the same time in the first incubation step. Secondary antibodies, which detect the primary antibodies, are added later. Therefore, the anti-FLAG antibody must be added with the other primary antibody.

    5.  **Conclusion:** Choice C correctly identifies both the timing (added with primary antibodies) and the specific, crucial reason for its inclusion (to verify GalNAc has not altered the surface expression of MUC1).
    """
    print(textwrap.dedent(explanation).strip())
    print("\nCorrect Answer Choice: C")

# Execute the function to provide the explanation.
explain_experimental_control()
