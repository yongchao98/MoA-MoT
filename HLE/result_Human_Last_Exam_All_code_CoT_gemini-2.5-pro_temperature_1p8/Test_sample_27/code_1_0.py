def solve_vj_puzzle():
    """
    This function determines and prints the significant contributors to the observation
    of dual light/alpha chains in single-cell V(D)J sequencing data for B and T cells.
    """

    # For B cells, significant contributors include technical artifacts (doublets, ambient RNA)
    # and key biological mechanisms related to light chain editing (a non-surface-expressed
    # chain or a remnant autoreactive chain). Having two fully functional BCRs is rare.
    b_cell_causes_indices = [1, 2, 4, 5]

    # For T cells, in addition to technical artifacts, the biology of TCR alpha chain
    # rearrangement makes dual expression very common. Therefore, cells with two functional
    # TCRs are a significant population, as are the related cases where one is non-functional
    # at the protein level or is autoreactive.
    t_cell_causes_indices = [1, 2, 3, 4, 5]

    # Format the numbers into the required string representation "(n,n,n)"
    # Using map(str,...) to convert each integer to a string before joining
    b_cell_causes_str = f"({','.join(map(str, b_cell_causes_indices))})"
    t_cell_causes_str = f"({','.join(map(str, t_cell_causes_indices))})"

    # Combine the two lists into the final output format: (list_B), (list_T)
    final_answer = f"{b_cell_causes_str}, {t_cell_causes_str}"

    print(final_answer)

solve_vj_puzzle()