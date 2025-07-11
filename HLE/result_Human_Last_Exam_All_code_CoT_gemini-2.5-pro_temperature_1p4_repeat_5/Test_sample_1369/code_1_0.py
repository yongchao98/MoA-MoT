import networkx as nx
import numpy as np

def analyze_graph():
    """
    This function builds a graph with 3 connected components and 2 cycles
    to demonstrate that this is a valid possibility under the problem's constraints,
    thus showing that we cannot conclude the graph has *exactly* two components.
    """
    # 1. Construct a counterexample graph with 3 components
    # Component 1: A triangle (3 vertices, 3 edges, 1 cycle)
    G1 = nx.cycle_graph(3)
    # Component 2: Another triangle (3 vertices, 3 edges, 1 cycle)
    G2 = nx.cycle_graph(3)
    # Component 3: A single edge (2 vertices, 1 edge, 0 cycles)
    G3 = nx.path_graph(2)

    # Create the full, disconnected graph
    G = nx.disjoint_union_all([G1, G2, G3])

    # 2. Verify its properties match the given information
    V = G.number_of_nodes()
    E = G.number_of_edges()
    k = nx.number_connected_components(G)

    # The cyclomatic number mu = |E| - |V| + k
    mu = E - V + k

    print("--- Analyzing a Counterexample Graph ---")
    print(f"Constructed a graph with {k} connected components.")
    print(f"Number of Vertices |V|: {V}")
    print(f"Number of Edges |E|: {E}")
    print(f"This demonstrates that a graph with k={k} components is possible.")
    
    print("\nVerifying the cyclomatic number condition:")
    print("The formula is: \u03BC = |E| - |V| + k")
    # Using the instruction to output each number in the final equation
    print(f"For our graph: {mu} = {E} - {V} + {k}")
    print(f"The result \u03BC={mu} matches the given condition (\u03BC=2).")

    # 3. Verify the eigenvalue information
    print("\nVerifying the eigenvalue condition:")
    L = nx.laplacian_matrix(G).toarray()
    eigenvalues = np.linalg.eigvalsh(L)
    # Count zero eigenvalues (with a small tolerance for floating point errors)
    zero_eigenvalues = np.sum(np.isclose(eigenvalues, 0))
    
    print(f"The number of zero eigenvalues is {zero_eigenvalues}, which equals k.")
    print("This is consistent with receiving the information that the first two eigenvalues are 0.")
    print("\n--- Conclusion ---")
    print("Since a graph with 3 components can satisfy all the given constraints,")
    print("we cannot conclude the graph *must* have exactly 2 components.")
    print("Therefore, none of the choices A, B, C, or D are guaranteed to be true.")

analyze_graph()