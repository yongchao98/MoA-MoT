import networkx as nx
import numpy as np
import scipy

def analyze_graph_property(G):
    """
    Analyzes a graph based on the property given in the problem.

    The property is null(B^T B) = lambda_n(G) / 2, which simplifies to
    c = lambda_n(G) / 2, where c is the cyclomatic number.
    """
    n = G.number_of_nodes()
    m = G.number_of_edges()
    k = nx.number_connected_components(G)

    print(f"Graph properties:")
    print(f"  Nodes (n): {n}")
    print(f"  Edges (m): {m}")
    print(f"  Connected components (k): {k}")
    
    # 1. Calculate the cyclomatic number, c = null(B^T B)
    c = m - n + k
    print(f"\n1. Calculate LHS of the equation:")
    print(f"   The cyclomatic number c = m - n + k")
    print(f"   c = {m} - {n} + {k} = {c}")
    print(f"   So, null(B^T B) = {c}")
    
    # 2. Calculate the largest Laplacian eigenvalue, lambda_n
    L = nx.laplacian_matrix(G).toarray()
    # Use eigvalsh for symmetric matrices
    eigenvalues = np.sort(scipy.linalg.eigvalsh(L))
    lambda_n = eigenvalues[-1]
    
    print(f"\n2. Calculate RHS of the equation:")
    print(f"   The largest Laplacian eigenvalue lambda_n = {lambda_n:.4f}")
    rhs = lambda_n / 2
    print(f"   lambda_n(G) / 2 = {lambda_n:.4f} / 2 = {rhs:.4f}")
    
    # 3. Check the identity
    print("\n3. Verify the identity null(B^T B) = lambda_n(G) / 2:")
    is_satisfied = np.isclose(c, rhs)
    print(f"   Is {c} == {rhs:.4f}? {'Yes' if is_satisfied else 'No'}")
    
    if not is_satisfied:
        print("This graph does not satisfy the initial condition.")
        return
        
    print("\n4. Check the answer choices for this graph:")
    # Choice B: k >= lambda_n / 2
    choice_b_holds = k >= rhs
    print(f"   B. k >= lambda_n/2  =>  {k} >= {rhs:.4f} ?  {choice_b_holds}")

    # Choice C: k == 1
    choice_c_holds = (k == 1)
    print(f"   C. is connected (k=1) ?  {choice_c_holds}")
    
    # Choice D: k == lambda_n / 2
    choice_d_holds = np.isclose(k, rhs)
    print(f"   D. k == lambda_n/2  =>  {k} == {rhs:.4f} ?  {choice_d_holds}")

# Create the graph G = C_4 U C_4
C4_1 = nx.cycle_graph(4)
C4_2 = nx.cycle_graph(4)
G = nx.disjoint_union(C4_1, C4_2)

# Run the analysis
analyze_graph_property(G)
