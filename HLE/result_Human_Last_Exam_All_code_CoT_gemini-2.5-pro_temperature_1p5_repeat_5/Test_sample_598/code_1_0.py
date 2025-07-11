import networkx as nx

def count_cycles_of_length_k(G, k):
    """
    Counts the number of simple cycles of length k in a graph G.
    
    This function works by finding all simple paths of length k-1 starting
    from each node, and then checking if an edge exists to close the cycle.
    To avoid overcounting, the final sum is divided by 2*k, as each cycle
    can be found starting from any of its k nodes, in either of 2 directions.
    
    Args:
        G (networkx.Graph): The graph to analyze.
        k (int): The length of the cycles to count.
        
    Returns:
        int: The total number of simple cycles of length k.
    """
    adj = {u: set(v for v in G.neighbors(u)) for u in G.nodes()}
    total_cycles = 0

    # Iterate over all nodes as potential starting points for a cycle
    for start_node in G.nodes():
        # q will store paths, processed layer by layer based on length.
        # Each item in q is a tuple: (path_list, path_set for fast lookups)
        q = [([start_node], {start_node})]
        
        # We build paths of length k-1 (which have k vertices)
        # We start with paths of length 0 and extend them k-1 times.
        for _ in range(k - 1):
            next_q = []
            for path, path_set in q:
                last_node = path[-1]
                for neighbor in adj[last_node]:
                    if neighbor not in path_set:
                        new_path = path + [neighbor]
                        new_path_set = path_set | {neighbor}
                        next_q.append((new_path, new_path_set))
            q = next_q
        
        # After the loop, q contains all simple paths of length k-1
        # starting at start_node. Now we check if any of them form a cycle.
        for path, path_set in q:
            last_node = path[-1]
            # Check for a closing edge from the last node back to the start node.
            # By construction, the start_node is not any of the intermediate
            # vertices in the path, so if this edge exists, it forms a k-cycle.
            if start_node in adj[last_node]:
                total_cycles += 1
                
    # Each cycle is counted 2k times (once per starting node in either direction)
    return total_cycles // (2 * k)

def main():
    """
    Constructs two non-isomorphic SRGs with parameters (16, 6, 2, 2),
    counts their 5-cycles, and prints the results.
    """
    print("This script demonstrates that two Strongly Regular Graphs (SRGs) with the same")
    print("parameters (n,d,lambda,mu) can have a different number of 5-cycles.\n")
    
    # Define the SRG parameters
    n, d, l, m = 16, 6, 2, 2
    k_cycle = 5
    
    print(f"We will examine two graphs from the class srg({n}, {d}, {l}, {m}).\n")

    # 1. Create the Rook's graph L(K_4,4)
    # This is equivalent to the Cartesian product of two K4 complete graphs.
    G1 = nx.cartesian_product(nx.complete_graph(4), nx.complete_graph(4))
    # NetworkX returns node labels as tuples, let's convert them to integers for clarity
    G1 = nx.convert_node_labels_to_integers(G1)
    
    # 2. Create the Shrikhande graph
    G2 = nx.shrikhande_graph()

    # Count 5-cycles in the first graph
    print("Calculating the number of 5-cycles for the first graph...")
    graph1_name = "Rook's graph on a 4x4 board (L(K_4,4))"
    c5_g1 = count_cycles_of_length_k(G1, k_cycle)
    
    print("Final Result for Graph 1:")
    print(f"Graph Name = {graph1_name}")
    print(f"Parameters (n,d,lambda,mu) = ({n}, {d}, {l}, {m})")
    print(f"Number of {k_cycle}-cycles = {c5_g1}\n")

    # Count 5-cycles in the second graph
    print("Calculating the number of 5-cycles for the second graph...")
    graph2_name = "Shrikhande graph"
    c5_g2 = count_cycles_of_length_k(G2, k_cycle)
    
    print("Final Result for Graph 2:")
    print(f"Graph Name = {graph2_name}")
    print(f"Parameters (n,d,lambda,mu) = ({n}, {d}, {l}, {m})")
    print(f"Number of {k_cycle}-cycles = {c5_g2}\n")

    # Conclusion
    if c5_g1 != c5_g2:
        print("Conclusion: The number of 5-cycles is different for these two graphs,")
        print(f"even though they share the same SRG parameters ({n},{d},{l},{m}).")
        print("Therefore, such a combination of parameters and graphs exists.")
    else:
        print("Conclusion: The number of 5-cycles is the same.")

if __name__ == '__main__':
    # The networkx library is required. Install with: pip install networkx
    main()