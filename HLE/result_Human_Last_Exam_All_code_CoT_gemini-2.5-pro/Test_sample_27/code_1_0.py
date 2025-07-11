import collections

def solve_vdj_observation():
    """
    This function determines and prints the significant contributors to observing
    dual light chains in B cells and dual alpha chains in T cells in a
    single-cell RNA-sequencing experiment.

    The reasoning is as follows:

    For B Cells (1 Heavy Chain, 2 Light Chains):
    - (1) Doublets: A significant technical artifact in droplet-based scRNA-seq.
    - (4) One transcript not on cell surface: A "failed" light chain that doesn't pair or traffic correctly can trigger another rearrangement. The original transcript can persist. This is a form of quality control and is significant.
    - (5) One transcript is autoreactive: This is the classic mechanism of receptor editing, a major biological process in B cell development to ensure self-tolerance. It is a very significant contributor.
    - (6) One transcript is not fully functional: Similar to (4), a non-functional light chain would be replaced. This is another facet of quality control.
    - Insignificant causes: (2) Measurement error is typically filtered, and (3) dual functional BCRs are very rare due to strong allelic exclusion.
    
    For T Cells (1 Beta Chain, 2 Alpha Chains):
    - (1) Doublets: A significant technical artifact.
    - (3) Both transcripts are functional and on the surface: This is a known biological phenomenon. Allelic exclusion at the TCR-alpha locus is inefficient, and up to 30% of T cells can express two functional TCRs.
    - (4) One transcript not on cell surface: During the process of generating a positively-selected TCR, multiple alpha chains may be transcribed. It's common for one to be functional while another is not expressed on the surface.
    - (6) One transcript is not fully functional: Similar to (4), this is a likely outcome of the continuous rearrangement process at the TCR-alpha locus.
    - Insignificant causes: (2) Measurement error is low, and (5) while a dual-alpha cell *could* be autoreactive, it's not the primary driving mechanism for dual expression as it is in B cells. T cells with high-affinity autoreactive TCRs are typically deleted.
    """

    # Define the significant causes for each cell type based on the reasoning above.
    # The numbers correspond to the mechanisms listed in the problem description.
    b_cell_causes = [1, 4, 5, 6]
    t_cell_causes = [1, 3, 4, 6]

    # Sort the lists for consistent ordering.
    b_cell_causes.sort()
    t_cell_causes.sort()
    
    # The problem asks us to output each number in the final equation.
    # We will format the output string as two comma-separated lists,
    # each enclosed by parentheses, separated by a comma and a space.
    b_cell_str = "(" + ",".join([str(num) for num in b_cell_causes]) + ")"
    t_cell_str = "(" + ",".join([str(num) for num in t_cell_causes]) + ")"
    
    final_answer = f"{b_cell_str}, {t_cell_str}"
    
    print(final_answer)

solve_vdj_observation()