import math
from itertools import combinations
from fractions import Fraction

def solve():
    """
    This function calculates the system-independent prefactor c_n for the
    virial expansion based on Ree-Hoover diagrams.
    """
    
    # The user can change this value to calculate c_n for a different n.
    # Note: The calculation time grows very rapidly with n. n=5 takes a few
    # seconds, n=6 can take several minutes.
    n = 4

    # --- Helper functions for graph analysis ---
    def _is_connected_util(num_vertices, adj_list):
        """Utility for checking connectivity from a starting node."""
        if num_vertices == 0:
            return True, 0
            
        start_node = 0
        
        q = [start_node]
        visited = {start_node}
        count = 1
        while q:
            u = q.pop(0)
            for v in adj_list[u]:
                if v not in visited:
                    visited.add(v)
                    q.append(v)
                    count += 1
        return count == num_vertices

    def is_connected(num_vertices, adj):
        """Checks if a graph is connected."""
        if num_vertices <= 1:
            return True
        
        # Check that there are no isolated vertices in the specified vertex set.
        # An articulation point check might create a disconnected graph of two non-trivial components.
        nodes_present = set()
        for i in range(num_vertices):
             if adj[i]:
                 nodes_present.add(i)
                 for neighbor in adj[i]:
                     nodes_present.add(neighbor)

        # if no edges, not connected for n > 1
        if not nodes_present:
            return num_vertices <= 1
        
        # Check connectivity of the non-isolated part of the graph
        # For our specific problem of articulation points, this means checking
        # if the remaining graph is connected.
        start_node = next(iter(nodes_present))
        q = [start_node]
        visited = {start_node}
        while q:
            u = q.pop(0)
            for v_neighbor in adj[u]:
                 if v_neighbor not in visited:
                    visited.add(v_neighbor)
                    q.append(v_neighbor)
        
        return visited == nodes_present


    def is_irreducible(num_vertices, edges):
        """
        Checks if a graph is irreducible.
        For n=2, irreducible means it's a single edge.
        For n>=3, irreducible means it is 2-connected.
        """
        if num_vertices <= 1:
            return False
        if num_vertices == 2:
            return len(edges) == 1

        # Build adjacency list
        adj = [set() for _ in range(num_vertices)]
        nodes_in_graph = set()
        for u, v in edges:
            adj[u].add(v)
            adj[v].add(u)
            nodes_in_graph.add(u)
            nodes_in_graph.add(v)
        
        # Subgraph must span all n vertices to be considered.
        if len(nodes_in_graph) != num_vertices:
            return False

        # 1. Check for basic connectivity
        adj_list_for_conn = [list(s) for s in adj]
        if not is_connected(num_vertices, adj_list_for_conn):
            return False

        # 2. Check for articulation points
        for i in range(num_vertices):
            # Create graph G' = G - {i}
            adj_prime = [[] for _ in range(num_vertices - 1)]
            
            def map_v(v):
                if v < i: return v
                if v > i: return v - 1
                return -1 

            for u_old in range(num_vertices):
                if u_old == i: continue
                u_new = map_v(u_old)
                for v_old in adj[u_old]:
                    if v_old == i: continue
                    v_new = map_v(v_old)
                    adj_prime[u_new].append(v_new)
            
            if not is_connected(num_vertices - 1, adj_prime):
                return False  # Found articulation point, so not 2-connected

        return True # No articulation points found, so 2-connected

    # --- Main calculation ---
    if n < 2:
        print("Error: n must be an integer >= 2.")
        return

    vertices = range(n)
    all_possible_edges = list(combinations(vertices, 2))
    num_edges_kn = len(all_possible_edges)

    c_kn = 0
    # Iterate over all 2^|E(K_n)| subgraphs
    for i in range(2**num_edges_kn):
        subgraph_edges = []
        for j in range(num_edges_kn):
            if (i >> j) & 1:
                subgraph_edges.append(all_possible_edges[j])

        if is_irreducible(n, subgraph_edges):
            num_edges_h = len(subgraph_edges)
            term = (-1)**(num_edges_kn - num_edges_h)
            c_kn += term

    n_minus_1 = n - 1
    n_factorial = math.factorial(n)
    
    cn_frac = -Fraction(n_minus_1 * c_kn, n_factorial)

    print(f"--- Calculation for c_{n} with n = {n} ---")
    print(f"The prefactor c_n is given by the formula: c_n = - (n-1)/n! * C(K_n)")
    print(f"For n = {n}, the components are:")
    print(f"  n-1 = {n_minus_1}")
    print(f"  n! = {n_factorial}")
    print(f"  The combinatorial factor C(K_{n}) = {c_kn}")
    print(f"Therefore, the final prefactor is:")
    print(f"  c_{n} = -({n_minus_1}) / {n_factorial} * ({c_kn}) = {cn_frac}")
    
    final_val_str = f"{cn_frac}"
    print(f"\n<<<{final_val_str}>>>")

solve()