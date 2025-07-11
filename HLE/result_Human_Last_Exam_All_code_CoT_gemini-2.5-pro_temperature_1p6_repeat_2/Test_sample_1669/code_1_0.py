import sys
import subprocess

# Ensure the pulp library is installed
try:
    import pulp
    import networkx as nx
except ImportError:
    print("pulp and networkx libraries are required. Please install them by running:")
    print("pip install pulp networkx")
    sys.exit(1)


def create_flower_snark_j5():
    """
    Constructs the Flower Snark J5 graph, which is a 3-regular graph with 20 vertices.
    This graph is known to be a "snark", meaning it's not 3-edge-colorable and lacks a 4-flow.
    """
    G = nx.Graph()
    # There are 4 groups of 5 vertices: u, v, w, x
    # We name them 'u0', 'u1', ... 'x4'
    for i in range(5):
        # The 'u' vertices form a 5-cycle
        G.add_edge(f'u{i}', f'u{(i + 1) % 5}')
        # The 'v' vertices form a 5-cycle with jumps of 2
        G.add_edge(f'v{i}', f'v{(i + 2) % 5}')
        # The 'w' vertices form a 5-cycle with jumps of -2 (or 3)
        G.add_edge(f'w{i}', f'w{(i + 3) % 5}')
        # Each 'x' vertex connects its corresponding u, v, and w vertices
        G.add_edge(f'x{i}', f'u{i}')
        G.add_edge(f'x{i}', f'v{i}')
        G.add_edge(f'x{i}', f'w{i}')
    return G

def find_k_vector(graph, k):
    """
    Tries to find a valid k-vector for the given graph using an Integer Linear Programming solver.
    A k-vector is an assignment of values from {+/-1, ..., +/-(k-1)} to edges
    such that the sum of values on edges incident to any vertex is 0.
    """
    print(f"Searching for a valid {k}-vector...")
    limit = k - 1
    if limit == 0:
        return None, f"k must be at least 2."

    # Use a canonical representation for edges (sorted tuple) to use as variable names
    edges = [tuple(sorted(edge)) for edge in graph.edges()]

    # 1. Create the ILP problem
    model = pulp.LpProblem(f"Valid_{k}-vector_search", pulp.LpMinimize) # Objective doesn't matter for feasibility

    # 2. Define variables
    # x_e is the integer value for each edge e, in the range [-limit, limit]
    x = pulp.LpVariable.dicts("x", edges, -limit, limit, pulp.LpInteger)

    # b_e is a binary variable to enforce x_e != 0
    b = pulp.LpVariable.dicts("b", edges, 0, 1, pulp.LpBinary)

    # 3. Add constraints
    # Constraint (a): For each vertex, the sum of values on incident edges must be 0
    for node in graph.nodes():
        incident_edges = [tuple(sorted(edge)) for edge in graph.edges(node)]
        model += pulp.lpSum(x[edge] for edge in incident_edges) == 0, f"flow_conservation_{node}"

    # Constraint (b): Enforce x_e is non-zero using the "big-M" method.
    # M=k is a sufficiently large constant.
    # if b[e]=0, then x[e] >= 1 (so x is positive)
    # if b[e]=1, then x[e] <= -1 (so x is negative)
    M = k
    for edge in edges:
        model += x[edge] >= 1 - M * b[edge], f"nonzero_lower_{edge}"
        model += x[edge] <= -1 + M * (1 - b[edge]), f"nonzero_upper_{edge}"

    # 4. Solve the model
    solver = pulp.PULP_CBC_CMD(msg=0)  # Use a quiet solver
    model.solve(solver)

    # 5. Return results
    status = pulp.LpStatus[model.status]
    if status == 'Optimal':
        # A solution was found
        solution = {edge: int(pulp.value(var)) for edge, var in x.items()}
        return solution, f"Success: A valid {k}-vector was found."
    else:
        # No solution exists
        return None, f"Failure: No valid {k}-vector exists for this graph."


def main():
    """Main function to run the analysis."""
    # We need to find the smallest k that works for ANY 20-vertex bridgeless 3-regular graph.
    # We test this on a "worst-case" example, the Flower Snark J5.
    j5_graph = create_flower_snark_j5()

    # Theory suggests k=4 will fail. Let's verify.
    solution_4, message_4 = find_k_vector(j5_graph, 4)
    print(message_4)
    
    if solution_4 is None:
        print("As predicted by theory, k=4 is not sufficient for all graphs.")
        # Theory suggests k=5 will succeed.
        solution_5, message_5 = find_k_vector(j5_graph, 5)
        print(message_5)

        if solution_5:
            smallest_k = 5
            print(f"\nThe smallest value of k is {smallest_k}.")

            # Display the "final equation" for one vertex as an example
            print("\nExample equation for vertex 'x0':")
            v = 'x0'
            incident_edges_vals = []
            for neighbor in j5_graph.neighbors(v):
                edge = tuple(sorted((v, neighbor)))
                val = solution_5[edge]
                incident_edges_vals.append(val)
            
            # Print the equation with signs
            equation_str = " + ".join(f"({v})" for v in incident_edges_vals).replace("+ (-", "- (")
            total = sum(incident_edges_vals)
            print(f"The values on the edges incident to '{v}' are: {incident_edges_vals[0]}, {incident_edges_vals[1]}, {incident_edges_vals[2]}")
            print(f"Equation: {equation_str} = {total}")

            # Final answer in the required format
            print("\n<<<5>>>")

if __name__ == '__main__':
    main()
