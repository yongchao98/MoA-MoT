def solve_vj_puzzle():
    """
    This function determines the significant contributors to the observation of dual
    light chains in B cells and dual alpha chains in T cells in a single-cell
    RNA-seq experiment.

    A mechanism is considered significant if it contributes more than 1% of the
    observed double light/alpha population.
    """

    # For B cells (one heavy chain, two light chains):
    # (1) Doublets: Significant technical artifact (~1-10% rate).
    # (2) Ambient RNA: Common source of noise, can be significant.
    # (3) True dual expression: Allelic exclusion is leaky for light chains; occurs in ~1-5% of B cells.
    # (4) One transcript non-surface: Non-productive rearrangements are common and can be transcribed.
    # (5) Receptor editing: A major biological driver for creating a second light chain to avoid autoreactivity.
    # (6) "Not fully functional" is too vague and overlaps with other, more precise causes.
    significant_mechanisms_b_cells = [1, 2, 3, 4, 5]

    # For T cells (one beta chain, two alpha chains):
    # (1) Doublets: Significant technical artifact, same as for B cells.
    # (2) Ambient RNA: Significant technical artifact, same as for B cells.
    # (3) True dual expression: Very common for TCR alpha chains (up to 30% of T cells) due to inefficient allelic exclusion.
    # (4) One transcript non-surface: Plausible result from the many rearrangement events.
    # (5) Autoreactive TCR: A known consequence of dual alpha expression, which complicates thymic selection. This describes a significant subset.
    # (6) "Not fully functional" is too vague.
    significant_mechanisms_t_cells = [1, 2, 3, 4, 5]

    # Format the lists into the required string format: (1,2,3), (4,5,6)
    # The map(str, ...) converts each number in the list to a string.
    # ", ".join(...) joins the string numbers with a comma and space.
    b_cell_string = f"({', '.join(map(str, significant_mechanisms_b_cells))})"
    t_cell_string = f"({', '.join(map(str, significant_mechanisms_t_cells))})"

    # Combine the two parts into the final answer string.
    final_answer = f"{b_cell_string}, {t_cell_string}"

    print(final_answer)


solve_vj_puzzle()
<<< (1, 2, 3, 4, 5), (1, 2, 3, 4, 5) >>>