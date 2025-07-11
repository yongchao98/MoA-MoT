import networkx as nx

def create_rooks_graph():
    """
    Constructs the 4x4 Rook's graph, which is an SRG with parameters (16, 6, 2, 2).
    Vertices are squares on a 4x4 board, edges connect squares in the same row or column.
    """
    n_size = 4
    G = nx.Graph()
    nodes = [(i, j) for i in range(n_size) for j in range(n_size)]
    G.add_nodes_from(nodes)
    for i in range(n_size):
        for j in range(n_size):
            for k in range(n_size):
                for l in range(n_size):
                    if (i, j) != (k, l):
                        if i == k or j == l:
                            G.add_edge((i, j), (k, l))
    return G

def create_shrikhande_graph():
    """
    Constructs the Shrikhande graph, also an SRG with parameters (16, 6, 2, 2).
    This is constructed as a Cayley graph on Z_4 x Z_4.
    """
    n_size = 4
    G = nx.Graph()
    nodes = [(i, j) for i in range(n_size) for j in range(n_size)]
    G.add_nodes_from(nodes)
    
    # Connection set for the Cayley graph
    # (0,+-1), (+-1,0), (1,1), (-1,-1)
    connection_set = {(0, 1), (0, 3), (1, 0), (3, 0), (1, 1), (3, 3)}

    for i in range(n_size):
        for j in range(n_size):
            u = (i, j)
            for dx, dy in connection_set:
                v = ((i + dx) % n_size, (j + dy) % n_size)
                G.add_edge(u, v)
    return G

def count_five_cycles(G):
    """
    Counts the number of simple 5-cycles in a graph.
    The method iterates through paths of length 4 and checks if the endpoints
    form a cycle. To avoid overcounting, the total sum is divided by 10
    (5 possible starting nodes, 2 directions for each cycle).
    """
    count = 0
    nodes = list(G.nodes())
    adj = G.adj
    
    for v1 in nodes:
        for v2 in adj[v1]:
            for v3 in adj[v2]:
                if v3 == v1:
                    continue
                for v4 in adj[v3]:
                    if v4 == v2 or v4 == v1:
                        continue
                    
                    # We have a path v1-v2-v3-v4 with distinct vertices.
                    # We need to close it to a 5-cycle with a vertex v5.
                    # v5 must be a common neighbor of v1 and v4.
                    
                    common_neighbors_v1_v4 = set(adj[v1]) & set(adj[v4])
                    for v5 in common_neighbors_v1_v4:
                        if v5 != v2 and v5 != v3:
                            count += 1
                            
    # Each 5-cycle (v1,v2,v3,v4,v5) is counted 2 * 5 = 10 times by this method.
    # v1 can be any of the 5 nodes, and the path can be traversed in 2 directions.
    return count // 10


def solve():
    """
    Main function to solve the problem.
    """
    params = {'n': 16, 'd': 6, 'lambda': 2, 'mu': 2}
    print(f"We investigate the class of strongly regular graphs srg(n,d,lambda,mu) for parameters:")
    print(f"n = {params['n']}, d = {params['d']}, lambda = {params['lambda']}, mu = {params['mu']}\n")
    
    # Create the two graphs
    rooks_graph = create_rooks_graph()
    shrikhande_graph = create_shrikhande_graph()
    
    # Count 5-cycles in each
    c5_rooks = count_five_cycles(rooks_graph)
    c5_shrikhande = count_five_cycles(shrikhande_graph)

    print("There are two non-isomorphic graphs for these parameters:")
    print("1. The Rook's graph (also known as the Latin square graph L_2(4))")
    print("2. The Shrikhande graph\n")
    
    print("We will now count the number of 5-cycles in each graph.")
    print(f"Number of 5-cycles in the Rook's graph: {c5_rooks}")
    print(f"Number of 5-cycles in the Shrikhande graph: {c5_shrikhande}\n")
    
    if c5_rooks != c5_shrikhande:
        print("As the counts are different, this provides a positive answer to the question.")
    else:
        print("The counts are the same.")

solve()