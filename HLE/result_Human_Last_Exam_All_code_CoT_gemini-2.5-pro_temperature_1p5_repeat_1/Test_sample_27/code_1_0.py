import collections

def solve_receptor_puzzle():
    """
    Solves the puzzle about multiple receptor chains in single-cell data.

    This function determines the significant contributors to the observation of B cells
    with one heavy chain and two light chains, and T cells with one beta chain and
    two alpha chains in a typical droplet-based scRNA-seq experiment.
    """

    # --- B Cell Analysis ---
    # In B cells, light chain allelic exclusion is imperfect or "leaky".
    #
    # Significant Contributors:
    # (1) Doublets: A major technical artifact where two cells are captured in one
    #     droplet. This is a very frequent source of multi-chain observations.
    # (3) True dual expression (both functional): It is well-documented that a
    #     fraction of B cells express two different functional light chains on their
    #     surface due to leaky allelic exclusion.
    # (5) Receptor editing: A key biological process where a B cell expressing an
    #     autoreactive receptor (Heavy + Light1) attempts to create a new, non-autoreactive
    #     light chain (Light2). The mRNA for the original Light1 chain can persist
    #     and be detected by scRNA-seq.
    # These three represent a major technical artifact and the two most prominent
    # biological mechanisms.
    b_cell_causes = collections.namedtuple('B_cell', ['technical_artifact', 'biological_dual_functional', 'biological_receptor_editing'])
    b_cell_indices = b_cell_causes(technical_artifact=1, biological_dual_functional=3, biological_receptor_editing=5)

    # --- T Cell Analysis ---
    # In T cells, there is strict allelic exclusion for the beta chain, but
    # there is NO allelic exclusion for the alpha chain.
    #
    # Significant Contributors:
    # (1) Doublets: As with B cells, this technical artifact is a significant source
    #     of apparent multi-alpha chain cells.
    # (3) True dual expression (both functional): The lack of alpha chain allelic
    #     exclusion is not a leak, but the normal mechanism. It commonly leads to
    #     T cells expressing two distinct, functional alpha chains, and thus two
    #     different TCRs. This is a major biological contributor, with up to 30% of
    #     peripheral T cells estimated to be dual-alpha.
    # (5) is not a mechanism for survival in T cells; autoreactivity leads to deletion.
    t_cell_causes = collections.namedtuple('T_cell', ['technical_artifact', 'biological_no_allelic_exclusion'])
    t_cell_indices = t_cell_causes(technical_artifact=1, biological_no_allelic_exclusion=3)

    # --- Formatting the output ---
    # The final output needs to show each number in the final list format:
    # (b_cell_list), (t_cell_list)

    # Convert tuples of indices to the required string format
    b_cell_output_list = [str(i) for i in b_cell_indices]
    t_cell_output_list = [str(i) for i in t_cell_indices]
    
    b_cell_output_str = f"({','.join(b_cell_output_list)})"
    t_cell_output_str = f"({','.join(t_cell_output_list)})"
    
    final_answer = f"{b_cell_output_str}, {t_cell_output_str}"

    print(final_answer)

solve_receptor_puzzle()
<<< (1,3,5), (1,3) >>>