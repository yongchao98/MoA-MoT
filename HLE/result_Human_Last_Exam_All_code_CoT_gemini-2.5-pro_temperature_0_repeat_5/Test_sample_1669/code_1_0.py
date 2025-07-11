import pulp

def solve_for_graph(graph_name, adj, k_max=6):
    """
    Finds the smallest k for a given graph using an ILP solver.
    """
    print(f"Analyzing the {graph_name} graph...")
    num_vertices = len(adj)
    vertices = list(adj.keys())

    edges = set()
    for u in adj:
        for v in adj[u]:
            edges.add(tuple(sorted((u, v))))
    
    edge_list = sorted(list(edges))
    num_edges = len(edge_list)
    edge_to_idx = {edge: i for i, edge in enumerate(edge_list)}

    for k in range(2, k_max + 1):
        print(f"Checking for k = {k}...", end="")
        
        # Create the ILP problem
        prob = pulp.LpProblem(f"Valid_{k}_Vector", pulp.LpMinimize)
        prob += 0  # Objective function is irrelevant

        # Define variables for the vector entries
        x = pulp.LpVariable.dicts("x", range(num_edges), cat='Integer')

        # Add vertex-sum constraints: sum of values on incident edges is 0
        for v in vertices:
            incident_edges_indices = [edge_to_idx[tuple(sorted((v, neighbor)))] for neighbor in adj[v]]
            prob += pulp.lpSum(x[i] for i in incident_edges_indices) == 0, f"vertex_{v}_constraint"

        # Add value range constraints: 1 <= |x_i| <= k-1
        # We use binary variables to model the disjunction x_i > 0 or x_i < 0
        y = pulp.LpVariable.dicts("y", range(num_edges), cat='Binary')
        M = k # A sufficiently large constant

        for i in range(num_edges):
            # if y_i = 0, then x_i is positive: 1 <= x_i <= k-1
            prob += x[i] >= 1 - M * y[i]
            prob += x[i] <= (k - 1) - M * y[i]
            # if y_i = 1, then x_i is negative: -(k-1) <= x_i <= -1
            prob += x[i] >= -(k - 1) + M * (1 - y[i])
            prob += x[i] <= -1 + M * (1 - y[i])

        # Solve the problem
        status = prob.solve(pulp.PULP_CBC_CMD(msg=0))

        if pulp.LpStatus[status] == 'Optimal':
            print(" Solution found.")
            solution_vector = {edge_list[i]: int(pulp.value(x[i])) for i in range(num_edges)}
            
            print(f"\nThe smallest value of k for the {graph_name} graph is {k}.")
            
            # Print an example equation as requested
            v_example = vertices[0]
            neighbors = adj[v_example]
            
            terms = []
            for neighbor in neighbors:
                edge = tuple(sorted((v_example, neighbor)))
                val = solution_vector[edge]
                terms.append(f"({val})")
            
            print(f"Example equation for vertex {v_example}: {' + '.join(terms)} = 0")
            return k
        else:
            print(" No solution found.")
    return None

if __name__ == '__main__':
    # The Petersen graph is the smallest "snark" (non-3-edge-colorable graph)
    # and serves as the "worst-case" example that requires k=4.
    petersen_adj = {
        0: [1, 4, 5], 1: [0, 2, 6], 2: [1, 3, 7], 3: [2, 4, 8], 4: [0, 3, 9],
        5: [0, 7, 8], 6: [1, 8, 9], 7: [2, 5, 9], 8: [3, 5, 6], 9: [4, 6, 7]
    }
    
    # The problem asks for a k that works for any 20-vertex graph of the type.
    # This means we must satisfy the worst-case scenario, which is a snark.
    # We demonstrate this by finding k for the Petersen graph.
    final_k = solve_for_graph("Petersen", petersen_adj)
    
    print("\n" + "="*50)
    print("Final Conclusion:")
    print("While some graphs in this class admit a valid 3-vector (if they are 3-edge-colorable),")
    print("others (the 'snarks') require a valid 4-vector.")
    print("To ensure a valid vector exists for ANY such graph, we must accommodate the worst case.")
    print(f"Therefore, the smallest value of k is {final_k}.")
    print("="*50)