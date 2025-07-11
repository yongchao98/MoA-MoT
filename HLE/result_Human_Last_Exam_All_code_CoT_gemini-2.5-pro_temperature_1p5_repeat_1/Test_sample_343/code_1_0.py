import networkx as nx
from itertools import combinations

def count_k_cliques(G, k):
    """
    Counts the number of k-cliques in a graph G.
    This is a brute-force algorithm.
    """
    count = 0
    nodes = list(G.nodes())
    for node_subset in combinations(nodes, k):
        subgraph = G.subgraph(node_subset)
        if subgraph.number_of_edges() == k * (k - 1) / 2:
            count += 1
    return count

def count_induced_k_matchings(G, k):
    """
    Counts the number of induced k-matchings in a graph G.
    This is a brute-force algorithm.
    """
    count = 0
    if k == 0:
        return 1
    
    edges = list(G.edges())
    for edge_subset in combinations(edges, k):
        vertices = set()
        is_matching = True
        for u, v in edge_subset:
            if u in vertices or v in vertices:
                is_matching = False
                break
            vertices.add(u)
            vertices.add(v)
        
        if is_matching:
            # Check if it's an induced matching
            subgraph = G.subgraph(vertices)
            if subgraph.number_of_edges() == k:
                count += 1
    return count

def count_induced_k_by_k_bicliques(G, k):
    """
    Counts the number of induced k-by-k bicliques in a graph G.
    This is a brute-force algorithm.
    """
    count = 0
    nodes = list(G.nodes())
    if 2 * k > len(nodes):
        return 0

    for node_subset in combinations(nodes, 2 * k):
        # Partition the 2k vertices into two sets A and B of size k
        # To avoid permutations and duplicates, we fix the smallest element to be in set A
        first_node = node_subset[0]
        remaining_nodes = node_subset[1:]
        
        for a_complement in combinations(remaining_nodes, k - 1):
            A = set([first_node] + list(a_complement))
            B = set(node_subset) - A
            
            # Check if A and B form an induced k,k-biclique
            # 1. No edges within A
            if G.subgraph(A).number_of_edges() != 0:
                continue
            
            # 2. No edges within B
            if G.subgraph(B).number_of_edges() != 0:
                continue

            # 3. All edges must exist between A and B
            is_complete = True
            for u in A:
                for v in B:
                    if not G.has_edge(u, v):
                        is_complete = False
                        break
                if not is_complete:
                    break
            
            if is_complete:
                count += 1
    return count


def main():
    """
    Main function to demonstrate the counting problem.
    """
    k = 3
    
    # Create a graph that is a disjoint union of:
    # 1. A 3-clique {0, 1, 2}
    # 2. An induced 3-matching on {3-8}
    # 3. An induced 3x3-biclique on {9-14}
    G = nx.Graph()

    # 1. Add a 3-clique
    G.add_edges_from([(0, 1), (0, 2), (1, 2)])
    
    # 2. Add an induced 3-matching
    G.add_edges_from([(3, 4), (5, 6), (7, 8)])

    # 3. Add an induced 3x3-biclique
    A = {9, 10, 11}
    B = {12, 13, 14}
    for u in A:
        for v in B:
            G.add_edge(u, v)

    print(f"Analyzing graph with {G.number_of_nodes()} nodes and {G.number_of_edges()} edges for k = {k}")

    num_cliques = count_k_cliques(G, k)
    num_induced_matchings = count_induced_k_matchings(G, k)
    num_induced_bicliques = count_induced_k_by_k_bicliques(G, k)
    
    total_count = num_cliques + num_induced_matchings + num_induced_bicliques
    
    print("\nCounting Results:")
    print(f"Number of {k}-cliques: {num_cliques}")
    print(f"Number of induced {k}-matchings: {num_induced_matchings}")
    print(f"Number of induced {k}-by-{k}-bicliques: {num_induced_bicliques}")
    
    print("\nFinal PCount equation:")
    print(f"{num_cliques} + {num_induced_matchings} + {num_induced_bicliques} = {total_count}")


if __name__ == "__main__":
    main()
