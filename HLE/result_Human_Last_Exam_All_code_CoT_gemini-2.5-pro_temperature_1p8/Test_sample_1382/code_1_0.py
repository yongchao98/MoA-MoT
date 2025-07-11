import networkx as nx
import numpy as np

def check_graph_property(G, name):
    """
    Checks a graph for the property specified in the problem.
    The property is: nullity(B^T B) = λ_n(G) / 2
    This translates to: m - n + c = λ_n / 2
    """
    n = G.number_of_nodes()
    m = G.number_of_edges()
    
    if n == 0:
        print(f"Graph: {name} is empty.")
        return

    c = nx.number_connected_components(G)
    
    # Cyclomatic number mu = m - n + c
    mu = m - n + c
    
    # Laplacian matrix and its eigenvalues
    L = nx.laplacian_matrix(G).toarray()
    eigenvalues = np.linalg.eigvalsh(L)
    lambda_n = eigenvalues[-1]

    # The two sides of the equation
    lhs = mu
    rhs = lambda_n / 2
    
    # Check the property given in the problem
    holds_prop = np.isclose(lhs, rhs)
    
    # Check the property from answer choice D
    holds_d = np.isclose(c, rhs)

    print(f"--- Analyzing graph: {name} ---")
    print(f"Nodes (n): {n}, Edges (m): {m}, Components (c): {c}")
    print(f"Left side of equation (μ = m - n + c): {lhs}")
    print(f"Right side of equation (λ_n / 2): {rhs:.4f}")
    
    print(f"Does the given property hold (μ == λ_n / 2)? {'Yes' if holds_prop else 'No'}")
    if holds_prop:
        print("  This suggests the graph is special.")
        
    print(f"Does answer D hold (c == λ_n / 2)? {'Yes' if holds_d else 'No'}")
    if np.isclose(m, n):
        print(f"Note: This graph has m=n, which implies μ=c.")
    else:
        print(f"Note: This graph has m!=n.")
    print("-" * 35 + "\n")


# --- Test Cases ---

# A graph where m != n
G1 = nx.complete_graph(4)
check_graph_property(G1, "K4 (Complete Graph, m > n)")

# A graph where m = n but the property doesn't hold
G2 = nx.cycle_graph(4)
check_graph_property(G2, "C4 (Cycle Graph, m = n)")

# A special graph (disjoint union of C3 and C4) where the property holds
C3 = nx.cycle_graph(3)
C4 = nx.cycle_graph(4)
G3 = nx.disjoint_union(C3, C4)
check_graph_property(G3, "C3 U C4 (Disconnected, m = n)")
