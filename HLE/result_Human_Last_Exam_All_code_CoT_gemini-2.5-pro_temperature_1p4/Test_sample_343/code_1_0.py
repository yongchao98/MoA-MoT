import itertools

def get_induced_subgraph(graph, nodes):
    """
    Computes the subgraph induced by a set of nodes.
    The graph is represented by an adjacency list.
    """
    subgraph = {u: [] for u in nodes}
    for u in nodes:
        for v in graph.get(u, []):
            if v in nodes:
                subgraph[u].append(v)
    return subgraph

def count_k_cliques(graph, k):
    """Counts the number of k-cliques in a graph."""
    if k <= 0:
        return 0
    nodes = list(graph.keys())
    n = len(nodes)
    if k > n:
        return 0
    
    count = 0
    for subset_nodes in itertools.combinations(nodes, k):
        is_clique = True
        for u, v in itertools.combinations(subset_nodes, 2):
            if v not in graph.get(u, []):
                is_clique = False
                break
        if is_clique:
            count += 1
    return count

def count_induced_k_matchings(graph, k):
    """Counts the number of induced k-matchings in a graph."""
    if k <= 0:
        return 0
    nodes = list(graph.keys())
    n = len(nodes)
    if 2 * k > n:
        return 0
        
    count = 0
    # Iterate over all subsets of 2k vertices
    for subset_nodes in itertools.combinations(nodes, 2 * k):
        subgraph = get_induced_subgraph(graph, set(subset_nodes))
        
        # An induced k-matching must have exactly k edges and all vertices must have degree 1
        num_edges = sum(len(neighbors) for neighbors in subgraph.values()) // 2
        
        if num_edges != k:
            continue
            
        is_induced_matching = True
        for node in subgraph:
            if len(subgraph[node]) != 1:
                is_induced_matching = False
                break
        
        if is_induced_matching:
            count += 1
    return count
    
def count_induced_k_by_k_bicliques(graph, k):
    """Counts the number of induced k-by-k bicliques in a graph."""
    if k <= 0:
        return 0
    nodes = list(graph.keys())
    n = len(nodes)
    if 2 * k > n:
        return 0

    count = 0
    # Iterate over all subsets of 2k vertices
    for subset_nodes in itertools.combinations(nodes, 2 * k):
        # Iterate over all ways to partition the 2k vertices into two sets of size k
        # The '/ 2' avoids double counting (A,B) and (B,A)
        for part_a_nodes in itertools.combinations(subset_nodes, k):
            part_a = set(part_a_nodes)
            part_b = set(subset_nodes) - part_a
            
            is_induced_biclique = True
            
            # Check for no edges within part_a
            for u, v in itertools.combinations(part_a, 2):
                if v in graph.get(u, []):
                    is_induced_biclique = False
                    break
            if not is_induced_biclique:
                continue

            # Check for no edges within part_b
            for u, v in itertools.combinations(part_b, 2):
                if v in graph.get(u, []):
                    is_induced_biclique = False
                    break
            if not is_induced_biclique:
                continue
            
            # Check for all edges between part_a and part_b
            for u in part_a:
                for v in part_b:
                    if v not in graph.get(u, []):
                        is_induced_biclique = False
                        break
                if not is_induced_biclique:
                    break
            
            if is_induced_biclique:
                count += 1
    
    # We divided by 2 when iterating partitions for a fixed set of 2k nodes.
    return count // 2


def main():
    # Example Graph G (a cycle on 6 vertices with one chord)
    # V = {0, 1, 2, 3, 4, 5}, E = {(0,1), (1,2), (2,3), (3,4), (4,5), (5,0), (0,3)}
    graph = {
        0: [1, 5, 3],
        1: [0, 2],
        2: [1, 3],
        3: [2, 4, 0],
        4: [3, 5],
        5: [4, 0]
    }
    
    # Let's take k = 3
    k = 3

    print(f"Solving PCount for the example graph with k = {k}\n")

    num_cliques = count_k_cliques(graph, k)
    num_matchings = count_induced_k_matchings(graph, k)
    num_bicliques = count_induced_k_by_k_bicliques(graph, k)
    
    total_count = num_cliques + num_matchings + num_bicliques
    
    print("The final equation is:")
    print(f"{num_cliques} (induced {k}-cliques) + "
          f"{num_matchings} (induced {k}-matchings) + "
          f"{num_bicliques} (induced {k}-by-{k}-bicliques) = "
          f"{total_count}")

if __name__ == '__main__':
    main()
