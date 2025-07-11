import math

def count_cliques_in_clique(n, k):
    """
    Calculates the number of k-cliques in a graph that is itself an n-clique.
    This is equivalent to the combination "n choose k".
    This function illustrates a component of a problem that is not FPT.
    Although C(n,k) is fast to compute for given n,k, finding all k-cliques 
    in a general graph takes O(N^k) time, where N is the number of vertices.
    This runtime is not FPT, supporting the claim that for patterns H with
    unbounded treewidth (like cliques where tw(K_k)=k-1), the problem is not FPT.
    """
    print(f"This helper function calculates C(n, k) for n={n}, k={k}.")
    if k < 0 or k > n:
        return 0
    if k == 0 or k == n:
        print("The final equation is trivial, C(n, k) = 1")
        print("The only number in the equation is: 1")
        return 1
    if k > n // 2:
        k = n - k
    
    # Build and print the full equation string as requested.
    num_terms = [str(i) for i in range(n, n - k, -1)]
    den_terms = [str(i) for i in range(k, 0, -1)]
    print(f"The equation is C({n}, {k}) = ({' * '.join(num_terms)}) / ({' * '.join(den_terms)})")

    # Print each number involved in the final equation.
    print("\nThe numbers involved in the numerator are:")
    for i in range(n, n-k, -1):
        print(i)
    
    print("\nThe numbers involved in the denominator are:")
    for i in range(k, 0, -1):
        print(i)
        
    # Calculate the result
    res = 1
    for i in range(k):
        res = res * (n - i) // (i + 1)
    
    return res

# --- Main block ---
# This is a demonstration based on the analysis.
n_nodes_in_clique = 20
k_subclique_size = 5
print(f"Counting {k_subclique_size}-cliques in a K_{n_nodes_in_clique} graph...")
print("-" * 30)
result = count_cliques_in_clique(n_nodes_in_clique, k_subclique_size)
print("-" * 30)
print(f"\nResult: There are {result} subgraphs isomorphic to K_{k_subclique_size} in a K_{n_nodes_in_clique}.")
