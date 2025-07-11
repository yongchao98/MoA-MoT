import collections
from itertools import combinations

def get_super_knight_graph(n, m):
    """
    Generates the super-knight graph for an n x m board.
    Returns an adjacency list representation of the graph.
    """
    adj = collections.defaultdict(list)
    moves = [
        (2, 3), (2, -3), (-2, 3), (-2, -3),
        (3, 2), (3, -2), (-3, 2), (-3, -2)
    ]
    for r in range(n):
        for c in range(m):
            u = r * m + c
            for dr, dc in moves:
                nr, nc = r + dr, c + dc
                if 0 <= nr < n and 0 <= nc < m:
                    v = nr * m + nc
                    adj[u].append(v)
    return adj

def find_k33_subgraph(adj):
    """
    Searches for a K_3,3 subgraph.
    Returns (U, V) if found, where U and V are lists of node indices.
    Otherwise, returns None.
    """
    nodes = sorted(adj.keys())
    if len(nodes) < 6:
        return None

    # Optimization: A vertex in a K_3,3 must have a degree of at least 3.
    # We will choose the first partition V from nodes with degree >= 3.
    candidate_v_nodes = [n for n in nodes if len(adj.get(n, [])) >= 3]
    if len(candidate_v_nodes) < 3:
        return None

    # Iterate through all combinations of 3 vertices for the V partition.
    for v_nodes in combinations(candidate_v_nodes, 3):
        # Find the common neighbors of the three vertices in V.
        # These common neighbors are candidates for the U partition.
        try:
            neighbors1 = set(adj[v_nodes[0]])
            neighbors2 = set(adj[v_nodes[1]])
            neighbors3 = set(adj[v_nodes[2]])
        except KeyError:
            continue
            
        common_neighbors = neighbors1.intersection(neighbors2, neighbors3)

        # If there are at least 3 common neighbors, we've found a K_3,3.
        if len(common_neighbors) >= 3:
            # U is a subset of size 3 of the common neighbors.
            u_nodes = list(combinations(list(common_neighbors), 3))[0]
            return (list(u_nodes), list(v_nodes))
            
    return None

def solve():
    """
    Finds the supremum of nm for which the super-knight graph is planar.
    """
    # Create a list of board sizes to check, sorted by area.
    sizes_to_check = []
    for n in range(4, 10):
        # Use m >= n to avoid checking both (4,5) and (5,4).
        for m in range(n, 10):
            sizes_to_check.append(((n, m), n * m))
    sizes_to_check.sort(key=lambda x: x[1])

    last_planar_area = 0

    print("Searching for the smallest non-planar super-knight graph...")
    print("-" * 60)

    for (n, m), area in sizes_to_check:
        print(f"Checking board size {n}x{m} (Area = {area})...", end="")
        adj = get_super_knight_graph(n, m)
        k33 = find_k33_subgraph(adj)

        if k33:
            print(" NON-PLANAR")
            print("-" * 60)
            print(f"Conclusion: The smallest non-planar graph is on a {n}x{m} board.")
            
            U_nodes, V_nodes = k33
            U_coords = sorted([(node // m, node % m) for node in U_nodes])
            V_coords = sorted([(node // m, node % m) for node in V_nodes])
            
            print(f"\nFound a K_3,3 subgraph as evidence:")
            print(f"Partition U: {U_coords}")
            print(f"Partition V: {V_coords}")
            
            supremum = last_planar_area
            print(f"\nThe largest area confirmed to be planar was {supremum}.")
            print(f"\nThe supremum of the value nm for which the graph G is planar is {supremum}.")
            return
        else:
            print(" Planar")
            last_planar_area = area
            
    print("Could not find a non-planar graph in the checked range.")

if __name__ == "__main__":
    solve()
<<<20>>>