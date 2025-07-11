import sys
from constraint import Problem, AllDifferentConstraint

# We need to increase Python's recursion limit for the constraint solver
sys.setrecursionlimit(2000)

def find_k_flow_for_petersen():
    """
    This function attempts to find a k-flow for the Petersen graph.
    The Petersen graph is the canonical example of a bridgeless 3-regular
    graph that does not have a 4-flow, thus requiring a 5-flow. This property
    is the reason why k=5 is the answer to the problem.

    We model this as a Constraint Satisfaction Problem.
    - Vertices: 0-9
    - Edges: 15. We define a variable for the flow on each edge.
    - Constraints: For each vertex, the sum of flows on incident edges must be zero (flow conservation).
    """

    # Petersen graph has 10 vertices and 15 edges.
    # We define the edges with a consistent orientation (u,v) where u < v.
    edges = [
        (0, 1), (0, 4), (0, 5), (1, 2), (1, 6), (2, 3), (2, 7),
        (3, 4), (3, 8), (4, 9), (5, 7), (5, 8), (6, 8), (6, 9), (7, 9)
    ]
    edge_vars = [f'e{i}' for i in range(len(edges))]

    # Adjacency list to easily define constraints
    adj = {i: [] for i in range(10)}
    for i, (u, v) in enumerate(edges):
        # The variable e_i represents flow from u to v.
        # This is an outgoing flow for u (-1) and incoming for v (+1).
        adj[u].append((f'e{i}', -1)) # Outgoing flow
        adj[v].append((f'e{i}', +1)) # Incoming flow

    for k in range(2, 6):
        print(f"--- Checking for a {k}-flow ---")
        problem = Problem()
        
        # A k-flow has non-zero integer values in {-(k-1), ..., -1, 1, ..., k-1}
        domain = list(range(-(k - 1), 0)) + list(range(1, k))
        problem.addVariables(edge_vars, domain)

        # Add flow conservation constraints for each vertex
        for vertex in range(10):
            # The lambda function captures the variables and coefficients for the sum
            # For each vertex, sum(coefficient * flow_variable) == 0
            variables = [var for var, c in adj[vertex]]
            coefficients = [c for var, c in adj[vertex]]
            problem.addConstraint(lambda *vals, c=coefficients: sum(c[i] * vals[i] for i in range(len(vals))) == 0, variables)
        
        # Find a solution
        solution = problem.getSolution()

        if solution:
            print(f"Success! Found a valid {k}-vector (a {k}-flow).")
            print("Flow values on each edge:")
            # Sort by edge variable name for consistent output
            sorted_vars = sorted(solution.keys())
            for var in sorted_vars:
                edge_index = int(var[1:])
                print(f"  Flow on edge {edges[edge_index]}: {solution[var]}")
            
            print("\nVerifying flow conservation at each vertex:")
            for vertex in range(10):
                terms = []
                total_flow = 0
                incident_flows = adj[vertex]
                for var, sign in incident_flows:
                    flow_value = solution[var] * sign
                    total_flow += flow_value
                    terms.append(f"({flow_value:+})")
                
                equation = f"V{vertex}: " + " + ".join(terms) + f" = {total_flow}"
                print(equation)
            
            print(f"\nThe smallest k for this graph is {k}.")
            print("\nBased on flow theory, since there exist 20-vertex graphs that inherit this")
            print("property from the Petersen graph, the smallest k for any such graph is 5.")
            return
        else:
            print(f"No solution found for k={k}. This is expected for k<5.")

if __name__ == '__main__':
    find_k_flow_for_petersen()

<<<5>>>