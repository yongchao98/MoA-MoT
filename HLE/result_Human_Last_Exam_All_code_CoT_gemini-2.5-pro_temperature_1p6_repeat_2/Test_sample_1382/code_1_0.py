import networkx as nx
import numpy as np

def get_nullity_unoriented_B(G):
    """
    Calculates the nullity of the unoriented incidence matrix B.
    nullity(B) = m - rank(B)
    """
    if G.number_of_edges() == 0:
        return 0
    # networkx gives B as a sparse matrix (n x m)
    B = nx.incidence_matrix(G, oriented=False).toarray()
    m = G.number_of_edges()
    rank_B = np.linalg.matrix_rank(B)
    return m - rank_B

def count_bipartite_non_tree_components(G):
    """
    Counts the number of connected components that are bipartite and not trees.
    """
    k = 0
    # Get component subgraphs
    component_subgraphs = [G.subgraph(c).copy() for c in nx.connected_components(G)]
    for comp in component_subgraphs:
        # A graph is a tree if it's connected and n = m + 1
        is_tree = nx.is_tree(comp)
        if nx.is_bipartite(comp) and not is_tree:
            k += 1
    return k

def main():
    """
    Main function to build the graph and check the properties.
    """
    # Create the graph G = C4 U C4 U K3
    G1 = nx.cycle_graph(4)
    G2 = nx.cycle_graph(4)
    G3 = nx.complete_graph(3)
    
    G = nx.disjoint_union(G1, G2)
    G = nx.disjoint_union(G, G3)
    
    # 1. Calculate the nullity of B^T B
    # This is equivalent to nullity(B)
    nullity_val = get_nullity_unoriented_B(G)
    
    # We can also calculate our property k directly for verification
    k = count_bipartite_non_tree_components(G)
    
    print(f"Graph G = C4 U C4 U K3 has n={G.number_of_nodes()} nodes and m={G.number_of_edges()} edges.")
    print("-" * 30)
    print(f"Verifying the interpretation of nullity(B^T B):")
    print(f"Directly computed nullity(B) = {nullity_val}")
    print(f"Count of bipartite, non-tree components (k) = {k}")
    print("-" * 30)

    # 2. Get the largest eigenvalue of the Laplacian
    L = nx.laplacian_matrix(G).toarray()
    eigenvalues = np.linalg.eigvalsh(L)
    lambda_n = np.max(eigenvalues)

    # 3. Check the condition given in the problem
    lhs = nullity_val
    rhs = lambda_n / 2
    
    print("Checking the problem's premise: null(B^T B) = lambda_n / 2")
    print(f"null(B^T B) = {lhs}")
    print(f"lambda_n / 2 = {lambda_n:.4f} / 2 = {rhs:.4f}")
    
    if np.isclose(lhs, rhs):
        print("The premise holds true for our example graph.")
    else:
        print("The premise does NOT hold true for our example graph.")
    print("-" * 30)
        
    # 4. Evaluate the answer choices
    c = nx.number_connected_components(G)
    
    print("Evaluating the answer choices for our graph:")
    # B. The graph has at least lambda_n(G)/2 connected components
    print("Choice B: c >= lambda_n / 2")
    print(f"Is {c} >= {rhs:.4f}?  --> {c >= rhs}")
    
    # C. The graph is connected
    print("\nChoice C: The graph is connected")
    print(f"Is {c} == 1?  --> {c == 1}")

    # D. The graph has exactly lambda_n(G)/2 connected components
    print("\nChoice D: c == lambda_n / 2")
    print(f"Is {c} == {rhs:.4f}?  --> {np.isclose(c, rhs)}")

if __name__ == '__main__':
    main()
