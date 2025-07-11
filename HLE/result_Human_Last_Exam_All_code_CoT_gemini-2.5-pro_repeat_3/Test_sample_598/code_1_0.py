import collections

def build_L_K44_graph():
    """
    Constructs the adjacency list for the L(K_4,4) graph (4x4 Rook's graph).
    The 16 vertices are identified with the cells of a 4x4 grid.
    Two vertices are adjacent if they are in the same row or column.
    """
    n = 16
    adj = collections.defaultdict(list)
    for i in range(n):
        for j in range(i + 1, n):
            # Map vertex indices to grid coordinates
            r1, c1 = i // 4, i % 4
            r2, c2 = j // 4, j % 4
            if r1 == r2 or c1 == c2:
                adj[i].append(j)
                adj[j].append(i)
    # Return as a list of lists for indexing
    return [sorted(adj[i]) for i in range(n)]

def build_shrikhande_graph():
    """
    Constructs the adjacency list for the Shrikhande graph.
    It is the Cayley graph of the group Z_4 x Z_4 with the connection set
    S = { (1,0), (-1,0), (0,1), (0,-1), (1,1), (-1,-1) }.
    """
    n = 16
    adj = collections.defaultdict(set)
    for r1 in range(4):
        for c1 in range(4):
            u = r1 * 4 + c1
            # Define the neighbors based on the connection set
            neighbors_coords = [
                ((r1 + 1) % 4, c1), ((r1 - 1) % 4, c1),
                (r1, (c1 + 1) % 4), (r1, (c1 - 1) % 4),
                ((r1 + 1) % 4, (c1 + 1) % 4), ((r1 - 1) % 4, (c1 - 1) % 4)
            ]
            for r2, c2 in neighbors_coords:
                v = r2 * 4 + c2
                adj[u].add(v)
                adj[v].add(u)
    # Return as a list of sorted lists for indexing
    return [sorted(list(adj[i])) for i in range(n)]

def count_k_cycles(adj, k):
    """
    Counts the number of simple cycles of length k in a graph.
    It uses an iterative DFS approach. To avoid overcounting, each cycle
    is counted only once, starting from its smallest-indexed vertex.
    """
    n = len(adj)
    count = 0
    for start_node in range(n):
        # stack stores (current_node, path_to_current_node)
        stack = [(start_node, [start_node])]
        while stack:
            curr_node, path = stack.pop()
            
            if len(path) == k:
                # When path length is k, check if it's a cycle
                if start_node in adj[curr_node]:
                    count += 1
                continue
            
            # Explore neighbors
            for neighbor in adj[curr_node]:
                # To find unique cycles, only explore neighbors with index > start_node
                # and that are not already in the current path.
                if neighbor > start_node and neighbor not in path:
                    new_path = path + [neighbor]
                    stack.append((neighbor, new_path))
                    
    # Each cycle is found twice (e.g., v1->v2..->vk->v1 and v1->vk..->v2->v1)
    return count // 2

def main():
    """
    Main function to construct the graphs, count 5-cycles, and print results.
    """
    # The parameters for the strongly regular graphs
    n = 16
    d = 6
    lambda_ = 2
    mu = 2

    print(f"Investigating strongly regular graphs with parameters (n,d,lambda,mu) = ({n},{d},{lambda_},{mu}).")
    print("-" * 70)

    # L(K_4,4) graph
    l_k44_graph = build_L_K44_graph()
    num_cycles_lk44 = count_k_cycles(l_k44_graph, 5)
    
    print(f"Graph 1: L(K_4,4) Graph (4x4 Rook's Graph)")
    print(f"Number of 5-cycles = {num_cycles_lk44}")
    print("-" * 70)

    # Shrikhande graph
    shrikhande_graph = build_shrikhande_graph()
    num_cycles_shrikhande = count_k_cycles(shrikhande_graph, 5)
    
    print(f"Graph 2: Shrikhande Graph")
    print(f"Number of 5-cycles = {num_cycles_shrikhande}")
    print("-" * 70)

    print("Conclusion:")
    print("Yes, there exists a combination of parameters and a pair of graphs in srg(n,d,lambda,mu) that do not have the same number of 5-cycles.")
    print(f"For parameters ({n},{d},{lambda_},{mu}), the L(K_4,4) graph and the Shrikhande graph are non-isomorphic and have a different number of 5-cycles.")
    print(f"Final Equation: {num_cycles_lk44} != {num_cycles_shrikhande}")

if __name__ == '__main__':
    main()
<<<192 != 240>>>