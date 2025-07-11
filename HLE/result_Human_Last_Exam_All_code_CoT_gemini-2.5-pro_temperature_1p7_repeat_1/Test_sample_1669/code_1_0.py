import pulp

def create_petersen_graph_edges():
    """Returns the edge list of the standard Petersen graph with vertices 0-9."""
    # Outer pentagon
    edges = [(i, (i + 1) % 5) for i in range(5)]
    # Inner star
    edges += [(i + 5, (i + 2) % 5 + 5) for i in range(5)]
    # Spokes
    edges += [(i, i + 5) for i in range(5)]
    return edges

def create_20_vertex_snark():
    """
    Constructs a 20-vertex snark by joining two Petersen graphs.
    This is done by taking two Petersen graphs, removing one edge from each,
    and then adding two new edges to connect the four resulting vertices of degree 2.
    The resulting graph is 3-regular, bridgeless, and not 3-edge-colorable.
    """
    # First Petersen graph (vertices 0-9)
    petersen1_edges = create_petersen_graph_edges()
    # Second Petersen graph (vertices 10-19, labels shifted by 10)
    petersen2_edges = [(u + 10, v + 10) for u, v in create_petersen_graph_edges()]

    # Remove edge (0, 1) from the first graph
    petersen1_edges.remove((0, 1))
    # Remove edge (10, 11) from the second graph
    petersen2_edges.remove((10, 11))
    
    # Vertices of the new graph
    vertices = list(range(20))
    # Edges of the new graph
    edges = petersen1_edges + petersen2_edges
    
    # Add new edges to connect the two modified graphs
    edges.append((0, 10))
    edges.append((1, 11))
    
    return vertices, edges

def find_k_flow(vertices, edges, k):
    """
    Tries to find a nowhere-zero k-flow for a given graph using an ILP solver.
    A k-flow has integer values on edges from the set {+/-1, ..., +/-(k-1)}.
    """
    if k <= 1:
        print("k must be > 1")
        return None

    # Create a directed graph for flow formulation
    oriented_edges = [(u, v) if u < v else (v, u) for u, v in edges]

    # Create the ILP problem
    prob = pulp.LpProblem("k-flow", pulp.LpMinimize) # Feasibility problem

    # ILP variables for flow on each edge
    # x_e can be negative, pulp handles this
    x = pulp.LpVariable.dicts("flow", oriented_edges, lowBound=-(k - 1), upBound=(k - 1), cat='Integer')

    # Binary variables to enforce non-zero flow
    b = pulp.LpVariable.dicts("is_neg", oriented_edges, cat='Binary')
    
    # Flow conservation constraints at each vertex
    for v in vertices:
        flow_sum = sum(x[e] for e in oriented_edges if e[1] == v) - \
                   sum(x[e] for e in oriented_edges if e[0] == v)
        prob += flow_sum == 0, f"flow_conservation_at_v{v}"

    # Non-zero constraints using the big-M method.
    # We want x_e in [-k+1, -1] U [1, k-1]
    M = k # A sufficiently large constant
    for e in oriented_edges:
        # if b[e] = 1, flow is negative. if b[e] = 0, flow is positive
        # Force x_e <= -1 if b[e] = 1
        prob += x[e] <= -1 + M * (1-b[e])
        # Force x_e >= 1 if b[e] = 0
        prob += x[e] >= 1 - M * b[e]

    # Solve the ILP
    print(f"Solving for a {k}-flow... (this may take a moment)")
    prob.solve(pulp.PULP_CBC_CMD(msg=0)) # Suppress solver output
    
    # Check the result
    if pulp.LpStatus[prob.status] == 'Optimal':
        print(f"\nSUCCESS: A valid {k}-vector (a nowhere-zero {k}-flow) was found.")
        flow_values = {e: int(x[e].varValue) for e in oriented_edges}
        print("A valid 5-vector assigns the following integer values to the edges of the graph:")
        # To avoid a very long output, we show the assignment for the first 10 edges
        # and summarize the set of values used.
        count = 0
        for edge, flow in flow_values.items():
            if count < 10:
                print(f"Edge {edge}: {flow}")
            count += 1
        if len(flow_values) > 10:
            print("...")
        
        unique_values = sorted(set(flow_values.values()))
        print(f"\nThe set of flow values used is: {unique_values}")
        print(f"All values are in the required range [{-k+1}, -1] U [1, {k-1}].")
        return flow_values
    else:
        print(f"\nFAILURE: No valid {k}-vector was found for this graph.")
        return None

if __name__ == '__main__':
    # 1. Define the graph
    V, E = create_20_vertex_snark()
    
    # 2. Set the value of k
    # We are testing k=5 because theory suggests it's the required value for this snark
    # and the minimal upper bound for the whole class of graphs.
    k_to_test = 5
    
    # 3. Find the k-vector
    find_k_flow(V, E, k_to_test)