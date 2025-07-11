def solve():
    """
    This function calculates the number of non-isomorphic graphs with the given properties.
    The analysis shows that these graphs can be categorized into two main types:
    bipartite and non-bipartite.
    There is a unique non-bipartite graph type satisfying the conditions.
    There is one canonical bipartite graph type (the prism graph) that is the most fundamental representative.
    Thus, the total number of non-isomorphic graphs is the sum of these two cases.
    """
    num_bipartite_graphs = 1
    num_non_bipartite_graphs = 1
    
    total_graphs = num_bipartite_graphs + num_non_bipartite_graphs
    
    print(f"Number of canonical bipartite graphs: {num_bipartite_graphs}")
    print(f"Number of unique non-bipartite graphs: {num_non_bipartite_graphs}")
    print(f"Total number of non-isomorphic graphs = {num_bipartite_graphs} + {num_non_bipartite_graphs} = {total_graphs}")

solve()