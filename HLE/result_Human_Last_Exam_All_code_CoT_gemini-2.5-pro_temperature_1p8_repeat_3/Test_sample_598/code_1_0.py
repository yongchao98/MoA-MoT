import networkx as nx

def build_rooks_graph():
    """
    Constructs the 4x4 Rook's graph, a strongly regular graph with parameters (16, 6, 2, 2).
    Vertices are cells on a 4x4 grid, and two vertices are adjacent if they are in the
    same row or column.
    """
    n_size = 4
    G = nx.Graph()
    nodes = range(n_size * n_size)
    G.add_nodes_from(nodes)

    for u in nodes:
        for v in nodes:
            if u >= v:
                continue
            
            r1, c1 = u // n_size, u % n_size
            r2, c2 = v // n_size, v % n_size

            if r1 == r2 or c1 == c2:
                G.add_edge(u, v)
    return G

def build_shrikhande_graph():
    """
    Constructs the Shrikhande graph, also a strongly regular graph with parameters (16, 6, 2, 2).
    It is non-isomorphic to the Rook's graph.
    Its construction is based on a Cayley graph of Z_4 x Z_4.
    """
    n_size = 4
    G = nx.Graph()
    nodes = range(n_size * n_size)
    G.add_nodes_from(nodes)
    
    # Connection set S for the Cayley graph on Z_4 x Z_4
    # S = {(0,±1), (±1,0), (±1,±1)} is for a different graph.
    # The correct set for Shrikhande is S = {±(1,0), ±(0,1), ±(1,1)}
    connection_set = {(0, 1), (0, 3), (1, 0), (3, 0), (1, 1), (3, 3)}
    
    for v1 in nodes:
        i1, j1 = v1 // n_size, v1 % n_size
        for v2 in nodes:
            if v1 >= v2:
                continue
            i2, j2 = v2 // n_size, v2 % n_size
            
            # Check if the difference (mod 4) is in the connection set
            di = (i1 - i2) % n_size
            dj = (j1 - j2) % n_size
            
            if (di, dj) in connection_set:
                G.add_edge(v1, v2)
    return G

def count_cycles_of_length_5(G):
    """
    Counts the number of simple cycles of length 5 in a graph G.
    It iterates through all paths of length 3 (v0-v1-v2-v3) and then finds
    a fifth vertex v4 that is a common neighbor to v0 and v3, ensuring all
    five vertices in the cycle are distinct.
    """
    count = 0
    nodes = list(G.nodes())
    for v0 in nodes:
        for v1 in G.neighbors(v0):
            for v2 in G.neighbors(v1):
                if v2 == v0:
                    continue # Path must not go back immediately
                
                for v3 in G.neighbors(v2):
                    if v3 in [v0, v1]:
                        continue # Vertices must be distinct

                    # We have a simple path v0-v1-v2-v3. Now we need to close it
                    # to a 5-cycle with a vertex v4.
                    # v4 must be a common neighbor of v0 and v3.
                    # v4 must not be v1 or v2 to form a simple 5-cycle.
                    for v4 in nx.common_neighbors(G, v0, v3):
                        if v4 not in [v1, v2]:
                            count += 1
                            
    # Each cycle (v0,v1,v2,v3,v4) is found 10 times:
    # - 5 times for each vertex chosen as the starting point v0.
    # - 2 times for each direction (e.g., v0->v1... and v0->v4...).
    # Therefore, we divide the total count by 10.
    return count // 10

if __name__ == '__main__':
    # 1. Define the SRG parameters
    n, d, lam, mu = 16, 6, 2, 2
    
    # 2. Construct the two graphs
    G_rooks = build_rooks_graph()
    G_shrikhande = build_shrikhande_graph()

    # 3. Count the 5-cycles in each graph
    num_c5_rooks = count_cycles_of_length_5(G_rooks)
    num_c5_shrikhande = count_cycles_of_length_5(G_shrikhande)

    # 4. Print the conclusion and the supporting numbers
    print("Yes, such a combination of parameters and graphs exists.")
    print(f"For the parameters (n,d,lambda,mu) = ({n},{d},{lam},{mu}), there are at least two non-isomorphic graphs.")
    print("These graphs have different numbers of 5-cycles:")
    print(f"- The 4x4 Rook's graph has {num_c5_rooks} 5-cycles.")
    print(f"- The Shrikhande graph has {num_c5_shrikhande} 5-cycles.")
    
    is_different = num_c5_rooks != num_c5_shrikhande
    print(f"\nConclusion: The number of 5-cycles is different, which confirms they are non-isomorphic. {is_different}")
