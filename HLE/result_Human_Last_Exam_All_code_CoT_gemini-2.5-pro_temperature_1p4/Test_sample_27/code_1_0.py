def solve_vdj_question():
    """
    This function determines the significant contributors to observing
    dual light/alpha chains in single-cell sequencing of B and T cells.

    The mechanisms are:
    (1) Cell doublets
    (2) Ambient RNA / measurement error
    (3) True dual expression, both functional on surface
    (4) True dual expression, one not on surface
    (5) True dual expression, one autoreactive
    (6) True dual expression, one not fully functional
    """

    # For B cells, technical artifacts (1,2) and biological mechanisms
    # related to receptor editing and leaky allelic exclusion (3,4,5,6)
    # are all significant contributors.
    b_cell_causes = (1, 2, 3, 4, 5, 6)

    # For T cells, technical artifacts (1,2) are significant. The biological
    # mechanism of dual alpha chain expression is very common due to inefficient
    # allelic exclusion, making all related scenarios (3,4,5,6) significant.
    t_cell_causes = (1, 2, 3, 4, 5, 6)

    # Format the answer as two comma-separated lists in parentheses.
    # The 'map(str, ...)' part ensures each number is converted to a string.
    # The '",".join(...)' part concatenates them with commas.
    b_cell_str = "(" + ",".join(map(str, b_cell_causes)) + ")"
    t_cell_str = "(" + ",".join(map(str, t_cell_causes)) + ")"

    # Combine the parts for the final answer string.
    final_answer = b_cell_str + ", " + t_cell_str

    print(final_answer)

solve_vdj_question()