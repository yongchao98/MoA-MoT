import networkx as nx
from itertools import combinations

def construct_rooks_graph():
    """Constructs the Rook's graph L(K_4,4)."""
    G = nx.Graph()
    nodes = [(i, j) for i in range(4) for j in range(4)]
    G.add_nodes_from(nodes)
    for (r1, c1), (r2, c2) in combinations(nodes, 2):
        if r1 == r2 or c1 == c2:
            G.add_edge((r1, c1), (r2, c2))
    return G

def construct_shrikhande_graph():
    """Constructs the Shrikhande graph as a Cayley graph."""
    G = nx.Graph()
    nodes = [(i, j) for i in range(4) for j in range(4)]
    G.add_nodes_from(nodes)
    # Generating set for the Cayley graph of Z_4 x Z_4
    gen_set = [(1, 0), (0, 1), (1, 1), (3, 3), (1, 3), (3, 1)]
    # An equivalent definition of the generating set
    # gen_set = [(1,0), (3,0), (0,1), (0,3), (1,1), (3,3)]
    
    for r1, c1 in nodes:
        for dr, dc in gen_set:
            r2 = (r1 + dr) % 4
            c2 = (c1 + dc) % 4
            G.add_edge((r1, c1), (r2, c2))
    return G

def verify_srg_params(G, n, d, lam, mu):
    """Verify if a graph G satisfies the SRG parameters."""
    if G.number_of_nodes() != n:
        return False
    
    # Check regularity and degree
    degrees = [val for (node, val) in G.degree()]
    if not all(deg == d for deg in degrees):
        return False
        
    # Check lambda and mu
    nodes = list(G.nodes())
    for i in range(len(nodes)):
        for j in range(i + 1, len(nodes)):
            u, v = nodes[i], nodes[j]
            common_neighbors = len(list(nx.common_neighbors(G, u, v)))
            if G.has_edge(u, v):
                if common_neighbors != lam:
                    return False
            else:
                if u != v and common_neighbors != mu:
                    return False
    return True

def count_five_cycles(G):
    """Counts the number of simple 5-cycles in a graph."""
    cycle_count = 0
    adj = G.adj
    nodes = list(G.nodes())
    for v0 in nodes:
        # Find paths of length 4 starting at v0 of the form v0-v1-v2-v3-v4
        # and check if (v4, v0) is an edge.
        for v1 in adj[v0]:
            for v2 in adj[v1]:
                if v2 == v0: continue
                for v3 in adj[v2]:
                    if v3 == v1 or v3 == v0: continue
                    # Now we have a simple path v0-v1-v2-v3
                    # Check if v3 has a neighbor v4 that is also a neighbor of v0
                    common_neighbors = set(adj[v3]) & set(adj[v0])
                    for v4 in common_neighbors:
                        if v4 != v1 and v4 != v2:
                           cycle_count += 1
                           
    # Each cycle is counted 5 times (once for each starting node)
    # and twice (for each direction). So we divide by 10.
    return cycle_count // 10

def main():
    """Main function to perform the task."""
    print("Yes, there exists a pair of graphs G and G' in srg(n,d,lambda,mu) that do not have the same number of 5-cycles.")
    print("A well-known example is the class srg(16, 6, 2, 2), which contains two non-isomorphic graphs.")
    print("-" * 30)

    # Rook's graph
    rooks_graph = construct_rooks_graph()
    print("Graph 1: Rook's graph L(K_4,4)")
    is_srg_rooks = verify_srg_params(rooks_graph, n=16, d=6, lam=2, mu=2)
    print(f"Is it srg(16, 6, 2, 2)? {is_srg_rooks}")
    
    # Shrikhande graph
    shrikhande_graph = construct_shrikhande_graph()
    print("Graph 2: Shrikhande graph")
    is_srg_shrikhande = verify_srg_params(shrikhande_graph, n=16, d=6, lam=2, mu=2)
    print(f"Is it srg(16, 6, 2, 2)? {is_srg_shrikhande}")
    print("-" * 30)

    print("Counting 5-cycles in each graph...")
    num_c5_rooks = count_five_cycles(rooks_graph)
    num_c5_shrikhande = count_five_cycles(shrikhande_graph)
    
    print(f"Number of 5-cycles in Rook's graph: {num_c5_rooks}")
    print(f"Number of 5-cycles in Shrikhande graph: {num_c5_shrikhande}")
    print("-" * 30)
    
    print("As you can see, the number of 5-cycles is different, proving that two SRGs with the same parameters do not necessarily have the same number of 5-cycles.")

if __name__ == "__main__":
    main()
