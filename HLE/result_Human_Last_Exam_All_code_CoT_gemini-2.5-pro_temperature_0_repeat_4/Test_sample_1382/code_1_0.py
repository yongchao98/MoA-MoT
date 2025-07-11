import numpy as np

def analyze_graph_condition():
    """
    Analyzes the problem by finding a counterexample.
    We test the complete graph K_4.
    """
    # Step 1: Define properties of K_4
    n = 4  # Number of nodes
    m = 6  # Number of edges
    c = 1  # Number of connected components
    # K_4 contains triangles (C_3), so it's not bipartite.
    # It has one connected component, which is not bipartite.
    c_b = 0 # Number of bipartite connected components

    print(f"Analyzing the complete graph K_4 with n={n}, m={m}, c={c}, c_b={c_b}\n")

    # Step 2: Calculate the left side of the condition: null(B^T B)
    # null(B^T B) = m - n + c_b
    null_B_T_B = m - n + c_b
    print(f"The term null(B^T B) is equal to m - n + c_b.")
    print(f"For K_4, this is {m} - {n} + {c_b} = {null_B_T_B}.")

    # Step 3: Calculate the right side of the condition: λ_n(L) / 2
    # Adjacency matrix for K_4
    A = np.ones((n, n)) - np.identity(n)
    # Degree matrix for K_4
    D = np.diag(np.sum(A, axis=1))
    # Laplacian matrix L = D - A
    L = D - A
    
    # Eigenvalues of the Laplacian
    eigenvalues = np.linalg.eigvalsh(L)
    lambda_n = np.max(eigenvalues)
    
    condition_right_side = lambda_n / 2
    
    print(f"\nThe Laplacian matrix for K_4 is:\n{L}")
    print(f"The eigenvalues of the Laplacian are: {np.round(eigenvalues, 2)}")
    print(f"The largest eigenvalue λ_n(L) is {lambda_n:.2f}.")
    print(f"The term λ_n(L) / 2 is {lambda_n:.2f} / 2 = {condition_right_side:.2f}.")

    # Step 4: Check if K_4 satisfies the condition
    print("\n--- Verifying the Condition ---")
    print(f"Is null(B^T B) = λ_n(L) / 2?")
    print(f"Is {null_B_T_B} = {condition_right_side}? {'Yes' if null_B_T_B == condition_right_side else 'No'}.")
    print("Since the condition holds, K_4 is a valid graph to test the answer choices.\n")

    # Step 5: Evaluate the answer choices for K_4
    print("--- Evaluating Answer Choices for K_4 ---")
    
    # Choice A
    print("A. If you drop λ_n(G)/2 edges, there will be at least two nodes with degree <= 1.")
    edges_to_drop = int(condition_right_side)
    print(f"   For K_4, λ_n(G)/2 = {edges_to_drop}. We need to drop {edges_to_drop} edges.")
    print("   K_4 has all degrees equal to 3. If we remove two disjoint edges, e.g., (0,1) and (2,3),")
    print("   the resulting graph is a 4-cycle (C_4), where all nodes have degree 2.")
    print("   In this case, no node has degree <= 1. So, statement A is FALSE.\n")

    # Choice B
    print("B. The graph has at least λ_n(G)/2 connected components.")
    print(f"   This means c >= λ_n(L) / 2.")
    print(f"   For K_4, is {c} >= {condition_right_side}? {c >= condition_right_side}.")
    print("   Statement B is FALSE.\n")

    # Choice C
    print("C. The graph is connected.")
    print(f"   For K_4, this is TRUE. However, we need a statement that is true for ALL graphs")
    print("   satisfying the condition. The empty graph on n>3 nodes also satisfies the condition")
    print("   (m=0, c=n, c_b=n, λ_n=0 => 0-n+n = 0/2), but it is not connected.")
    print("   Therefore, statement C is FALSE in general.\n")

    # Choice D
    print("D. The graph has exactly λ_n(G)/2 connected components.")
    print(f"   This means c = λ_n(L) / 2.")
    print(f"   For K_4, is {c} = {condition_right_side}? {c == condition_right_side}.")
    print("   Statement D is FALSE.\n")

    print("--- Conclusion ---")
    print("Since we found a valid graph (K_4) for which choices A, B, and D are false,")
    print("and another valid graph (the empty graph) for which C is false,")
    print("none of the statements A, B, C, or D are necessary consequences of the condition.")

analyze_graph_condition()