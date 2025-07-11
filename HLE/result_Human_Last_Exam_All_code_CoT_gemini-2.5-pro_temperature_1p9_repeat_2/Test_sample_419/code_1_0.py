def solve_biology_question():
    """
    Analyzes an experimental biology question and provides the correct answer with an explanation.
    """
    explanation = """
The experiment's goal is to demonstrate that an antibody specifically binds to the TN antigen on MUC1, using a competitive inhibition assay with free N-acetylgalactosamine (GalNAc). A key challenge with using a high concentration of an inhibitor like 500 mM GalNAc is the potential for non-specific, or 'off-target', effects on the cells. A primary concern is that such a high concentration could cause cellular stress, leading to the internalization or shedding of surface proteins, including the MUC1 construct.

If MUC1 is lost from the cell surface due to this stress, the anti-MUC1 antibody would fail to bind, mimicking the effect of specific competitive inhibition and leading to an incorrect conclusion.

To control for this, it is essential to measure the level of MUC1 protein on the cell surface independently of its glycosylation. The MUC1 construct carries a FLAG tag for this exact purpose. An anti-FLAG antibody will bind to this tag, providing a reliable measure of MUC1 surface expression.

This control has two requirements:
1.  **Purpose:** It must verify that MUC1 surface expression is unchanged by the 500 mM GalNAc treatment.
2.  **Timing:** The anti-FLAG antibody is a 'primary' antibody because it binds directly to its target on the MUC1 protein. Therefore, it should be added during the primary antibody incubation step, along with the main antibody of interest (the anti-TN MUC1 antibody).

Based on this analysis:
-   Choice C correctly identifies that the anti-FLAG antibody should be added with the primary antibodies.
-   Choice C provides the most precise and critical reason: to verify that the high concentration of GalNAc has not altered the surface expression of MUC1.
"""
    print(explanation)
    final_answer = "C"
    # No calculation needed, but the prompt says to explain the thinking. The text above does that.

solve_biology_question()