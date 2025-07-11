import networkx as nx

def calculate_poset_global_dimension(poset_relations):
    """
    Calculates the global dimension of the incidence algebra of a finite poset.

    Args:
        poset_relations (list of tuples): A list of pairs (a, b) representing the
                                         relation a < b in the poset.
    """
    if not poset_relations:
        print("Global dimension for an empty or discrete poset is 0.")
        return

    # Create a directed graph from the poset relations (edges from smaller to larger)
    G = nx.DiGraph()
    G.add_edges_from(poset_relations)
    
    # Ensure the graph is a valid poset (acyclic)
    if not nx.is_directed_acyclic_graph(G):
        raise ValueError("The provided relations do not form a valid poset (contains cycles).")

    # Get a topological sort of the nodes to process them from minimal to maximal
    try:
        nodes = list(nx.topological_sort(G))
    except nx.NetworkXUnfeasible:
        # This case is caught by the is_directed_acyclic_graph check, but good practice
        raise ValueError("The provided relations do not form a valid poset.")

    # Memoization cache for projective dimensions of simple modules
    pd = {}

    print("Calculating projective dimensions for simple modules S_i:")
    for i in nodes:
        # Find all predecessors of node i
        predecessors = list(G.predecessors(i))

        if not predecessors:
            # i is a minimal element
            pd[i] = 0
            print(f"  pd(S_{i}) = 0 (since {i} is a minimal element)")
        else:
            # i is not minimal. Find the maximal predecessors.
            # An element j is a maximal predecessor if no other predecessor k satisfies j < k.
            maximal_predecessors = []
            for j in predecessors:
                is_maximal = True
                for k in predecessors:
                    if j != k and nx.has_path(G, j, k):
                        is_maximal = False
                        break
                if is_maximal:
                    maximal_predecessors.append(j)
            
            # The projective dimension of S_i is 1 + max(pd of simples for maximal predecessors)
            max_pd_of_predecessors = max(pd[j] for j in maximal_predecessors)
            pd[i] = 1 + max_pd_of_predecessors
            print(f"  pd(S_{i}) = 1 + max{{pd(S_{j}) for j in {maximal_predecessors}}} = 1 + {max_pd_of_predecessors} = {pd[i]}")
    
    # Global dimension is the maximum of all these projective dimensions
    if not pd:
        global_dim = 0
    else:
        global_dim = max(pd.values())
        
    pd_values = list(pd.values())
    print("\n---")
    print("Final Calculation:")
    print(f"Global Dimension = max{{pd(S_i) for all i}}")
    print(f"Global Dimension = max({pd_values})")
    print(f"Global Dimension = {global_dim}")

# Example 1: A chain C_4 = {1 < 2 < 3 < 4}. (Tame, finite type)
# Its global dimension should be 3.
print("--- Example 1: Chain poset C_4 = {1 < 2 < 3 < 4} ---")
chain_poset = [(1, 2), (2, 3), (3, 4)]
calculate_poset_global_dimension(chain_poset)

# Example 2: The 4-subspace poset. (Tame, non-finite type)
# Its global dimension should be 1.
print("\n--- Example 2: The 4-subspace poset (1,1,1,1) ---")
four_subspace_poset = [('a', 'm'), ('b', 'm'), ('c', 'm'), ('d', 'm')]
calculate_poset_global_dimension(four_subspace_poset)

# Example 3: Poset J from {1<3, 2<3}, whose incidence algebra is not hereditary. (Rep-finite)
# Its global dimension should be 2.
print("\n--- Example 3: The poset J where {1 < 3, 2 < 3} ---")
v_poset = [(1, 3), (2, 3)]
calculate_poset_global_dimension(v_poset)
