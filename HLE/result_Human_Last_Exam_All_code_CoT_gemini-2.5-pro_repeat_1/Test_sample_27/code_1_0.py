def solve_receptor_puzzle():
    """
    This function determines the significant contributors to observing
    dual light/alpha chains in single-cell V(D)J sequencing data.

    Based on the analysis of both technical artifacts and biological phenomena:
    - B-Cells: All 6 mechanisms are considered significant.
    - T-Cells: All 6 mechanisms are considered significant.
    """

    # For B cells, doublets (1), ambient RNA (2), allelic inclusion (3),
    # non-productive transcripts (4), receptor editing (5), and ongoing
    # selection (6) are all significant contributors.
    b_cell_contributors = (1, 2, 3, 4, 5, 6)

    # For T cells, doublets (1), ambient RNA (2), the lack of alpha-chain
    # allelic exclusion (3), non-productive transcripts (4), escape from
    # negative selection (5), and ongoing thymic selection (6) are all
    # significant contributors.
    t_cell_contributors = (1, 2, 3, 4, 5, 6)

    # Format the numbers into the required string format "(n,n,...)"
    b_cell_str = f"({','.join(map(str, b_cell_contributors))})"
    t_cell_str = f"({','.join(map(str, t_cell_contributors))})"

    # Combine the two parts into the final answer string
    final_answer = f"{b_cell_str}, {t_cell_str}"

    print(final_answer)

solve_receptor_puzzle()
<<< (1,2,3,4,5,6), (1,2,3,4,5,6) >>>