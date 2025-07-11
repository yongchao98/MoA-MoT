def solve_vj_puzzle():
    """
    This function determines and prints the significant contributors to the observation
    of dual light/alpha chains in B/T cells in single-cell RNA-seq data.
    """

    # Significant contributors for B cells with 1 heavy chain and 2 light chains
    # 1: Doublets (technical artifact)
    # 2: Ambient RNA (technical artifact)
    # 4: Dual mRNA, one not on surface (biological: leaky exclusion)
    # 5: Dual mRNA, one autoreactive (biological: receptor editing)
    b_cell_contributors = (1, 2, 4, 5)

    # Significant contributors for T cells with 1 beta chain and 2 alpha chains
    # 1: Doublets (technical artifact)
    # 2: Ambient RNA (technical artifact)
    # 3: Dual surface expression (biological: no TCRa allelic exclusion)
    # 4: Dual mRNA, one not on surface (biological: no TCRa allelic exclusion)
    t_cell_contributors = (1, 2, 3, 4)

    # Format the numbers into strings with comma separation
    b_cell_str = ", ".join(map(str, b_cell_contributors))
    t_cell_str = ", ".join(map(str, t_cell_contributors))

    # Construct the final answer string in the specified format
    final_answer = f"({b_cell_str}), ({t_cell_str})"

    print(final_answer)

solve_vj_puzzle()