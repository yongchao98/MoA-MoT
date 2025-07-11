import math

def solve_researcher_problem():
    """
    Solves the combinatorics problem about researchers at a conference.
    """
    # Step 1: Define the initial parameters from the problem statement.
    N = 42  # Total number of researchers
    d = 24  # Number of collaborators for each researcher (degree of each vertex)
    N0 = 2027 # Number of tables with 0 collaborations

    # Step 2: Set up the first equation by counting the total number of possible tables.
    # Total tables = C(N, 3) = N0 + N1 + N2 + N3
    total_tables = math.comb(N, 3)
    # This gives us: N1 + N2 + N3 = total_tables - N0
    sum_N1_N2_N3 = total_tables - N0

    print(f"Total number of researchers: {N}")
    print(f"Collaborators per researcher: {d}")
    print(f"Number of tables with no collaborations (N0): {N0}")
    print(f"Total possible ways to form a table of 3: C(42, 3) = {total_tables}")
    print("\n--- Equation 1: Based on Total Tables ---")
    print(f"N1 + N2 + N3 = {total_tables} - {N0}")
    print(f"N1 + N2 + N3 = {sum_N1_N2_N3}\n")

    # Step 3: Set up the second equation by counting the total number of edges in all tables.
    # Total edges in the graph = (N * d) / 2
    total_edges_in_graph = (N * d) // 2
    # Each edge is part of (N - 2) tables.
    # Total edges in all tables = total_edges_in_graph * (N - 2)
    sum_edges_in_all_tables = total_edges_in_graph * (N - 2)
    # This sum also equals 1*N1 + 2*N2 + 3*N3
    
    print("--- Equation 2: Based on Total Edges ---")
    print(f"Total collaboration pairs (edges) in the graph = (42 * 24) / 2 = {total_edges_in_graph}")
    print(f"Each edge is in (42 - 2) = {N-2} tables.")
    print(f"Sum of edges over all possible tables = {total_edges_in_graph} * {N-2} = {sum_edges_in_all_tables}")
    print(f"N1 + 2*N2 + 3*N3 = {sum_edges_in_all_tables}\n")

    # Step 4: Set up the third equation by counting total "V-shapes" (paths of length 2).
    # V-shapes centered at one vertex = C(d, 2)
    v_shapes_per_vertex = math.comb(d, 2)
    # Total V-shapes = N * C(d, 2)
    total_v_shapes = N * v_shapes_per_vertex
    # This sum also equals 1*N2 + 3*N3
    
    print("--- Equation 3: Based on Total 'V-shapes' ---")
    print(f"V-shapes centered at one researcher = C(24, 2) = {v_shapes_per_vertex}")
    print(f"Total V-shapes in the graph = 42 * {v_shapes_per_vertex} = {total_v_shapes}")
    print(f"N2 + 3*N3 = {total_v_shapes}\n")

    # Step 5: Solve the system of linear equations for N3.
    # System:
    # (1) N1 + N2 + N3 = sum_N1_N2_N3
    # (2) N1 + 2*N2 + 3*N3 = sum_edges_in_all_tables
    # (3)      N2 + 3*N3 = total_v_shapes

    # Subtract (1) from (2) to get a new equation (4) with N2 and N3.
    # (4) N2 + 2*N3 = sum_edges_in_all_tables - sum_N1_N2_N3
    sum_N2_2N3 = sum_edges_in_all_tables - sum_N1_N2_N3
    
    print("--- Solving the System ---")
    print("By subtracting Equation 1 from Equation 2, we get:")
    print(f"N2 + 2*N3 = {sum_edges_in_all_tables} - {sum_N1_N2_N3} = {sum_N2_2N3}\n")

    # Subtract (4) from (3) to solve for N3.
    # N3 = total_v_shapes - sum_N2_2N3
    N3 = total_v_shapes - sum_N2_2N3
    
    print("By subtracting this new equation from Equation 3, we can find N3:")
    print(f"N3 = ({total_v_shapes}) - ({sum_N2_2N3})")
    print(f"The number of constellations where all three researchers have authored with each other is {N3}.")

solve_researcher_problem()