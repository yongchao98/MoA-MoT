def solve_sc_rsv_puzzle():
    """
    This function determines the significant contributors to the observation of
    dual light/alpha chains in single-cell B/T cell sequencing data and
    prints the result in the specified format.
    """

    # For B cells, significant contributors include technical errors (doublets, ambient RNA)
    # and biological phenomena (leaky allelic exclusion, receptor editing to remove autoreactivity).
    b_cell_causes = [1, 2, 3, 4, 5]

    # For T cells, significant contributors include technical errors (doublets, ambient RNA)
    # and biological phenomena (common dual alpha chain expression, non-productive transcripts).
    # Autoreactive cells are clonally deleted, so this is not a significant cause.
    t_cell_causes = [1, 2, 3, 4, 6]

    # Format the lists into the required string output format.
    # The format is (list1), (list2)
    b_cell_str = ",".join(map(str, b_cell_causes))
    t_cell_str = ",".join(map(str, t_cell_causes))

    final_answer_string = f"({b_cell_str}), ({t_cell_str})"

    print(final_answer_string)

solve_sc_rsv_puzzle()