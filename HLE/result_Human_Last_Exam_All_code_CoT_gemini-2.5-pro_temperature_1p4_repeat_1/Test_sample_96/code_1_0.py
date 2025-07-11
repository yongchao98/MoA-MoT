import collections
from itertools import permutations

def solve_e8_torsion_problem():
    """
    This script solves the problem by counting pairs of disconnected subgraphs (A1, A4)
    within the E8 Dynkin diagram.
    """
    
    # The E8 Dynkin diagram represented by an adjacency list.
    # Nodes are labeled 1 to 8.
    #       2
    #       |
    # 1--3--4--5--6--7--8
    e8_adj = {
        1: {3}, 2: {4}, 3: {1, 4}, 4: {2, 3, 5},
        5: {4, 6}, 6: {5, 7}, 7: {6, 8}, 8: {7}
    }
    nodes = set(e8_adj.keys())

    # Step 1: Find all subgraphs of type A4 (paths on 4 vertices).
    # We iterate through all permutations of 4 nodes and check if they form a path.
    a4_subgraphs = set()
    for p in permutations(nodes, 4):
        # A permutation (n1, n2, n3, n4) is a path if n1-n2, n2-n3, n3-n4 are edges.
        if p[1] in e8_adj[p[0]] and p[2] in e8_adj[p[1]] and p[3] in e8_adj[p[2]]:
            # Store the set of nodes as a frozenset to handle duplicates.
            a4_subgraphs.add(frozenset(p))

    # Step 2: For each node `i` (representing an A1 subgroup), count how many
    # A4 subgraphs `J` are disconnected from it.
    counts_per_node = collections.OrderedDict()
    total_count = 0
    
    for i in sorted(list(nodes)):
        count_for_i = 0
        # A subgraph is disconnected from node `i` if it contains no nodes
        # adjacent to `i`.
        neighbors_of_i = e8_adj.get(i, set())
        
        for j_nodes in a4_subgraphs:
            # The A4 subgraph `J` cannot contain node `i` itself.
            if i in j_nodes:
                continue
            
            # Check if any node in `J` is a neighbor of `i`.
            if j_nodes.isdisjoint(neighbors_of_i):
                count_for_i += 1
        
        counts_per_node[i] = count_for_i
        total_count += count_for_i

    # Step 3: Print the results as requested.
    # The total number of elements is the sum of the counts for each node.
    equation_parts = [str(v) for v in counts_per_node.values()]
    print(f"The number of elements can be calculated by summing the valid configurations for each node representing the A1 subgroup:")
    print(f"{' + '.join(equation_parts)} = {total_count}")


solve_e8_torsion_problem()