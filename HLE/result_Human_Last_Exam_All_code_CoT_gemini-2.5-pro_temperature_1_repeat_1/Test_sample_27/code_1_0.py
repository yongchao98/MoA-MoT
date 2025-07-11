def solve_vj_puzzle():
    """
    This function determines and prints the significant contributors to the observation
    of single cells with one heavy/beta chain and two light/alpha chains for B and T cells.
    """

    # For B cells, significant contributors include technical artifacts (doublets, ambient RNA)
    # and biological phenomena (receptor editing, leaky allelic exclusion).
    # (1) Doublets: Two cells in one droplet.
    # (2) Ambient RNA: Contamination from lysed cells.
    # (3) Dual Surface Expression: Result of leaky exclusion.
    # (4) One Non-Surface Chain: Result of leaky exclusion or editing.
    # (5) Receptor Editing: Primary mechanism to remove autoreactivity.
    b_cell_causes = (1, 2, 3, 4, 5)

    # For T cells, receptor editing (5) does not occur. The primary biological cause
    # is the naturally inefficient allelic exclusion of the TCR alpha chain.
    # (1) Doublets: Two cells in one droplet.
    # (2) Ambient RNA: Contamination from lysed cells.
    # (3) Dual Surface Expression: Common result of leaky alpha chain exclusion.
    # (4) One Non-Surface Chain: Plausible result of leaky alpha chain exclusion.
    t_cell_causes = (1, 2, 3, 4)

    # Format the output string as per the requirements:
    # Two comma-separated lists, each in parentheses, with B cells first.
    b_cell_str = ", ".join(map(str, b_cell_causes))
    t_cell_str = ", ".join(map(str, t_cell_causes))

    final_answer = f"({b_cell_str}), ({t_cell_str})"
    print(final_answer)

solve_vj_puzzle()