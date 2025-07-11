def solve_vj_puzzle():
    """
    This function determines and prints the significant contributors to observing
    B cells with two light chains and T cells with two alpha chains in scRNA-seq.
    """

    # Significant contributors for B cells
    # (1) Doublets (technical artifact)
    # (2) Ambient RNA (technical artifact)
    # (4) One productive, one non-productive transcript (biological reality)
    # (5) One autoreactive transcript being replaced via receptor editing (biological reality)
    b_cell_causes = [1, 2, 4, 5]

    # Significant contributors for T cells
    # (1) Doublets (technical artifact)
    # (2) Ambient RNA (technical artifact)
    # (3) Two functional surface-expressed transcripts due to leaky allelic exclusion (biological reality)
    # (4) One productive, one non-productive transcript (biological reality)
    t_cell_causes = [1, 2, 3, 4]

    # Format the output as two comma-separated lists in parentheses
    b_cell_str = "(" + ", ".join(map(str, b_cell_causes)) + ")"
    t_cell_str = "(" + ", ".join(map(str, t_cell_causes)) + ")"

    final_answer = b_cell_str + ", " + t_cell_str
    
    print(final_answer)

solve_vj_puzzle()