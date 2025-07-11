def solve_vj_puzzle():
    """
    This function determines and prints the significant contributors to the observation
    of B cells with one heavy and two light chains, and T cells with one beta and
    two alpha chains in a single-cell RNA-seq experiment.

    The final answer is formatted as two comma-separated lists, each enclosed
    by parentheses, and separated by a comma and a space.
    """

    # For B cells, significant contributors are:
    # (1) Technical error: doublets (two cells in one droplet).
    # (2) Technical error: ambient RNA contamination / barcode swapping.
    # (4) Biological reason: As a result of receptor editing, one of the light chain transcripts
    #     is for a receptor that is being replaced and will not be expressed on the surface.
    # (5) Biological reason: Receptor editing is initiated when the first receptor is autoreactive.
    #     This is the primary driver for generating a second light chain.
    b_cell_contributors = (1, 2, 4, 5)

    # For T cells, significant contributors are:
    # (1) Technical error: doublets (two cells in one droplet).
    # (2) Technical error: ambient RNA contamination / barcode swapping.
    # (3) Biological reason: Allelic exclusion for the TCR alpha chain is inefficient,
    #     leading to a large population of T cells with two functional alpha chains on the surface.
    t_cell_contributors = (1, 2, 3)

    # Format the output string as required
    b_cell_str = f"({', '.join(map(str, sorted(b_cell_contributors)))})"
    t_cell_str = f"({', '.join(map(str, sorted(t_cell_contributors)))})"
    
    # Print the final answer
    final_answer = f"{b_cell_str}, {t_cell_str}"
    print(final_answer)

solve_vj_puzzle()
# The final output will be printed to the console.
# We expect the output to be: (1, 2, 4, 5), (1, 2, 3)