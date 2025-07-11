# The user needs to install the following libraries for the code to run:
# pip install networkx pulp

import networkx as nx
import pulp

def find_smallest_k_for_graph(G):
    """
    Finds the smallest k for which a graph G has a nowhere-zero k-flow
    by solving an integer linear program.
    """
    nodes = list(G.nodes())
    edges = list(G.edges())
    edge_to_idx = {edge: i for i, edge in enumerate(edges)}

    # Create an arbitrary orientation for edges (e.g., from smaller to larger node index)
    # and an incidence mapping for vertices.
    incidence = {node: [] for node in nodes}
    for i, edge in enumerate(edges):
        u, v = edge
        source, target = (u, v) if u < v else (v, u)
        incidence[source].append((i, -1))  # Sign -1 for outgoing flow
        incidence[target].append((i, +1))  # Sign +1 for incoming flow

    # Loop through k to find the smallest value that admits a flow.
    # We test up to k=6 based on Seymour's 6-flow theorem.
    for k in range(2, 7):
        print(f"--- Checking for a valid {k}-vector (k={k}) ---")

        # Create the ILP problem. We want a feasible solution, so the objective is trivial.
        prob = pulp.LpProblem(f"Nowhere_Zero_{k}-Flow_Problem", pulp.LpMinimize)
        prob += 0, "Arbitrary Objective"

        # Define integer variables for the flow on each edge.
        # Flow values can be in {-(k-1), ..., -1, 1, ..., k-1}.
        flow_vars = pulp.LpVariable.dicts("flow", range(len(edges)), cat='Integer',
                                          lowBound=-(k - 1), upBound=(k - 1))
        
        # Add binary variables to enforce the "nowhere-zero" constraint.
        binary_vars = pulp.LpVariable.dicts("binary", range(len(edges)), cat='Binary')

        # Add constraints to the problem.
        # 1. Flow conservation at each node (sum of flows is zero).
        for node in nodes:
            prob += pulp.lpSum(sign * flow_vars[edge_idx] for edge_idx, sign in incidence[node]) == 0, f"Flow_Conservation_at_Node_{node}"

        # 2. Nowhere-zero flow constraint for each edge using the Big-M method.
        # This forces flow_vars[i] to be either >= 1 or <= -1.
        M = k  # A sufficiently large constant
        for i in range(len(edges)):
            prob += flow_vars[i] >= 1 - M * binary_vars[i], f"Non_Zero_Positive_Edge_{i}"
            prob += flow_vars[i] <= -1 + M * (1 - binary_vars[i]), f"Non_Zero_Negative_Edge_{i}"

        # Solve the ILP. Suppress verbose output from the solver.
        solver = pulp.PULP_CBC_CMD(msg=False)
        prob.solve(solver)

        # Check the result.
        status = pulp.LpStatus[prob.status]
        print(f"Result: {status}")

        if status == 'Optimal':
            print(f"\nSuccess! A valid {k}-vector exists.")
            
            # Print the flow conservation equation for the first vertex as an example.
            node_to_check = nodes[0]
            print(f"Example Equation (flow conservation at node {node_to_check}):")
            
            equation_numbers = []
            for edge_idx, sign in incidence[node_to_check]:
                # The term in the sum is sign * flow_value
                flow_value = flow_vars[edge_idx].varValue
                term = sign * flow_value
                equation_numbers.append(str(int(round(term))))

            equation_str = " + ".join(f"({n})" for n in equation_numbers)
            print(f"{equation_str} = 0")
            
            return k
        print("-" * 35)

    return None

if __name__ == '__main__':
    # We use the Flower Snark J5, a known 20-vertex, 3-regular, bridgeless graph
    # that requires a 5-flow.
    G = nx.flower_snark()

    print("Analyzing graph properties to confirm it matches the description:")
    print(f"Graph name: Flower Snark J5")
    print(f"Number of vertices: {G.number_of_nodes()}")
    print(f"Number of edges: {G.number_of_edges()}")
    is_3_regular = all(d == 3 for n, d in G.degree())
    print(f"Is 3-regular? {is_3_regular}")
    is_bridgeless = nx.is_bridgeless(G)
    print(f"Is bridgeless? {is_bridgeless}\n")

    # Run the main function to find k
    smallest_k = find_smallest_k_for_graph(G)

    if smallest_k is not None:
        print(f"\nConclusion: The smallest value of k such that the graph admits a valid k-vector is {smallest_k}.")

<<<5>>>