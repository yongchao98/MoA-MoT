import collections

def find_significant_contributors():
    """
    This function identifies the significant contributors to the observation of dual
    light/alpha chains in B/T cells in a scRNA-seq experiment.

    The function will:
    1. Define the lists of significant cause-identifiers for B cells and T cells based on immunological
       and technical principles.
    2. Format the two lists into the specified output string: (list_B), (list_T)
    3. Print the final result.
    """

    # For B cells, expressing one heavy chain and two light chains:
    # (1) Doublets and (2) Ambient RNA are significant technical artifacts.
    # (3) True dual functional surface receptors are very rare and selected against. Not significant.
    # (4, 5, 6) Leaky allelic exclusion and receptor editing are significant biological
    # mechanisms that result in two light chain mRNAs.
    b_cell_causes = [1, 2, 4, 5, 6]
    b_cell_causes.sort()

    # For T cells, expressing one beta chain and two alpha chains:
    # (1) Doublets and (2) Ambient RNA are significant technical artifacts.
    # (3, 4, 5, 6) The continuous and non-exclusive rearrangement of the TCR alpha locus
    # makes all these biological scenarios significant contributors. Dual TCR-alpha
    # expression is a known, common phenomenon.
    t_cell_causes = [1, 2, 3, 4, 5, 6]
    t_cell_causes.sort()

    # Format the lists into strings as per the requested format.
    # For example, [1, 2, 4] becomes "(1,2,4)"
    b_cell_str = "(" + ",".join(map(str, b_cell_causes)) + ")"
    t_cell_str = "(" + ",".join(map(str, t_cell_causes)) + ")"

    # Combine the two strings for the final output.
    final_answer_string = f"{b_cell_str}, {t_cell_str}"

    print(final_answer_string)

find_significant_contributors()