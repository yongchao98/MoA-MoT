import itertools

def get_connected_subsets(nodes, edges):
    """Find all non-empty connected subsets of a graph."""
    adj = {u: set() for u in nodes}
    for u, v in edges:
        adj[u].add(v)
        adj[v].add(u)

    all_connected_subsets = set()
    for i in range(1, len(nodes) + 1):
        for subset_nodes in itertools.combinations(nodes, i):
            subset_nodes = frozenset(subset_nodes)
            
            # Check for connectivity using BFS/DFS
            q = [next(iter(subset_nodes))]
            visited = {q[0]}
            head = 0
            while head < len(q):
                u = q[head]
                head += 1
                for v in adj[u]:
                    if v in subset_nodes and v not in visited:
                        visited.add(v)
                        q.append(v)
            
            if visited == subset_nodes:
                all_connected_subsets.add(subset_nodes)
                
    return all_connected_subsets

def apply_automorphism(subset, perm):
    """Apply a graph automorphism (permutation) to a subset."""
    return frozenset(perm[node] for node in subset)

def count_orbits(subsets, automorphisms):
    """Count the number of orbits of subsets under a group of automorphisms."""
    unclassified = set(subsets)
    num_orbits = 0
    while unclassified:
        num_orbits += 1
        representative = unclassified.pop()
        
        orbit = set()
        for perm in automorphisms:
            orbit.add(apply_automorphism(representative, perm))
            
        unclassified -= orbit
    return num_orbits

def main():
    """
    Calculates the number of distinct types of connected subgraphs for
    a path graph (analogue of interval [0,1]) and a cycle graph (analogue of circle S^1).
    """
    # 1. Path Graph P_4 (analogue for [0,1])
    p4_nodes = frozenset([0, 1, 2, 3])
    p4_edges = [(0, 1), (1, 2), (2, 3)]
    p4_subsets = get_connected_subsets(p4_nodes, p4_edges)
    # Automorphisms of P4: identity and reflection
    p4_autos = [
        {0: 0, 1: 1, 2: 2, 3: 3},  # identity
        {0: 3, 1: 2, 2: 1, 3: 0},  # reflection
    ]
    p4_orbits = count_orbits(p4_subsets, p4_autos)
    
    print("Analysis for Path Graph P4 (analogue for interval [0,1]):")
    print(f"Total number of connected subgraphs: {len(p4_subsets)}")
    print(f"Number of distinct types (orbits): {p4_orbits}\n")

    # 2. Cycle Graph C_3 (analogue for S^1)
    c3_nodes = frozenset([0, 1, 2])
    c3_edges = [(0, 1), (1, 2), (2, 0)]
    c3_subsets = get_connected_subsets(c3_nodes, c3_edges)
    # Automorphisms of C3: Dihedral group D3 (rotations and reflections)
    c3_autos = [
        {0: 0, 1: 1, 2: 2},  # id
        {0: 1, 1: 2, 2: 0},  # rot 1
        {0: 2, 1: 0, 2: 1},  # rot 2
        {0: 0, 1: 2, 2: 1},  # reflect across 0
        {0: 2, 1: 1, 2: 0},  # reflect across 1
        {0: 1, 1: 0, 2: 2},  # reflect across 2
    ]
    c3_orbits = count_orbits(c3_subsets, c3_autos)

    print("Analysis for Cycle Graph C3 (analogue for circle S^1):")
    print(f"Total number of connected subgraphs: {len(c3_subsets)}")
    print(f"Number of distinct types (orbits): {c3_orbits}\n")

    final_answer = c3_orbits
    print(f"The smallest number of topologically distinct compactifications is: {final_answer}")


if __name__ == "__main__":
    main()
