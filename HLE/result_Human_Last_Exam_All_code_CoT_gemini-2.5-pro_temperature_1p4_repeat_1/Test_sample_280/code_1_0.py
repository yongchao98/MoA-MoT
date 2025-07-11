import math

def solve_researcher_problem():
    """
    Calculates the number of table constellations where all three researchers
    have authored papers with each other.
    """
    # Step 1: Define the initial parameters from the problem description.
    N = 42  # Total number of researchers
    d = 24  # Number of co-authors for each researcher (degree of each vertex)
    E0_count = 2027  # Number of tables with 0 co-author pairs (independent sets of size 3)

    # Step 2: Calculate the values needed for the system of equations.

    # Total number of possible tables (triplets), C(N, 3)
    total_tables = math.comb(N, 3)
    
    # Equation 1: Sum of all table categories
    # |E_1| + |E_2| + |E_3| = total_tables - |E_0|
    sum_E1_E2_E3 = total_tables - E0_count
    
    # Equation 2: Total number of "wedges"
    # |E_2| + 3 * |E_3| = N * C(d, 2)
    num_wedges = N * math.comb(d, 2)

    # Equation 3: Total number of (edge, table) pairs
    # |E_1| + 2*|E_2| + 3*|E_3| = |Edges| * (N - 2)
    num_edges = (N * d) // 2
    edge_table_pairs = num_edges * (N - 2)

    # Step 3: Solve the system of equations for |E_3|.
    # We have:
    # (A) |E_1| + |E_2| +   |E_3| = sum_E1_E2_E3
    # (B)         |E_2| + 3*|E_3| = num_wedges
    # (C) |E_1| + 2*|E_2| + 3*|E_3| = edge_table_pairs

    # Subtracting (A) from (C) gives:
    # (D) |E_2| + 2*|E_3| = edge_table_pairs - sum_E1_E2_E3
    sum_E2_2E3 = edge_table_pairs - sum_E1_E2_E3

    # Subtracting (D) from (B) gives:
    # |E_3| = num_wedges - sum_E2_2E3
    E3_count = num_wedges - sum_E2_2E3

    # Step 4: Print the final calculation and result.
    print("The number of tables where all three researchers are co-authors, |E_3|, is calculated as follows:")
    print(f"Total tables with at least one co-author pair (|E_1|+|E_2|+|E_3|) = C({N}, 3) - {E0_count} = {total_tables} - {E0_count} = {sum_E1_E2_E3}")
    print(f"Total 'wedges' (|E_2| + 3*|E_3|) = {N} * C({d}, 2) = {num_wedges}")
    print(f"Total edge-table pairs (|E_1| + 2*|E_2| + 3*|E_3|) = ({N}*{d}/2) * ({N}-2) = {edge_table_pairs}")
    print("\nSolving the system of equations leads to the final calculation for |E_3|:")
    print(f"|E_3| = (Total wedges) - (Total edge-table pairs - Total tables with at least one co-author pair)")
    print(f"|E_3| = {num_wedges} - ({edge_table_pairs} - {sum_E1_E2_E3})")
    print(f"|E_3| = {num_wedges} - {sum_E2_2E3}")
    print(f"|E_3| = {E3_count}")
    
    print(f"\nThus, the number of table constellations where all three researchers have authored with each other is {E3_count}.")

solve_researcher_problem()