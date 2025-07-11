def solve_task():
    """
    This function determines the significant contributors for observing
    cells with one heavy/beta chain and two light/alpha chains in
    single-cell RNA sequencing data.
    """

    # For B cells, significant contributors include technical artifacts (doublets, ambient RNA)
    # and biological realities (failed allelic exclusion, receptor editing).
    b_cell_causes = [1, 2, 3, 4, 5]

    # For T cells, the same technical artifacts apply. Biologically, dual alpha chain
    # expression is a very common and expected feature of T cell development.
    t_cell_causes = [1, 2, 3, 4, 5]

    # Format the output string as requested, e.g., (1,2,3), (4,5,6)
    b_cell_str = "(" + ",".join(map(str, b_cell_causes)) + ")"
    t_cell_str = "(" + ",".join(map(str, t_cell_causes)) + ")"

    final_answer = b_cell_str + ", " + t_cell_str
    print(final_answer)

solve_task()