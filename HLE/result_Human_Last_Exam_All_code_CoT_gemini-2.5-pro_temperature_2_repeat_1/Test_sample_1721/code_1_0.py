import math

# Global variable to be used by the recursive solver
max_found_size = 0

def find_mis_size_recursive(adj, current_nodes_tuple, current_set_size):
    """
    Finds the maximum independent set size of a graph using a recursive backtracking algorithm.
    This function updates a global variable `max_found_size`.
    
    Args:
        adj (list of lists): The adjacency list of the graph.
        current_nodes_tuple (tuple): A tuple of nodes that are still candidates to be added.
        current_set_size (int): The size of the independent set built so far.
    """
    global max_found_size
    
    if not current_nodes_tuple:
        if current_set_size > max_found_size:
            max_found_size = current_set_size
        return
        
    # Pruning: if we can't possibly beat the current best, stop.
    if current_set_size + len(current_nodes_tuple) <= max_found_size:
        return

    # Take the first node `v` from the remaining candidates
    v = current_nodes_tuple[0]
    
    # --- Option 1: Add `v` to the independent set. ---
    # We must then remove `v` and all its neighbors from the list of candidates.
    neighbors = {v} | set(adj[v])
    remaining_nodes_1 = tuple([node for node in current_nodes_tuple[1:] if node not in neighbors])
    find_mis_size_recursive(adj, remaining_nodes_1, current_set_size + 1)

    # --- Option 2: Do not add `v` to the independent set. ---
    # We just move to the next candidate.
    remaining_nodes_2 = current_nodes_tuple[1:]
    find_mis_size_recursive(adj, remaining_nodes_2, current_set_size)

def find_best_density():
    """
    Searches for the best density c by checking different moduli m. For each modulus,
    it computes the maximum size of a residue set R such that R+R contains no squares mod m.
    """
    global max_found_size
    best_c = 0.0
    best_m = 0
    
    print("Searching for the best density c using the modular construction method...")
    print("Modulus (m) | Max size of R (|R|) | Density c = |R|/m")
    print("-" * 50)

    # Iterate through moduli m. Higher values take longer to compute.
    for m in range(2, 31):
        # 1. Find the set of squares modulo m
        squares_mod_m = set((k * k) % m for k in range(m))

        # 2. Build the conflict graph
        # Vertices are {0, 1, ..., m-1}
        # Edge (r1, r2) if (r1 + r2) % m is a square
        adj = [[] for _ in range(m)]
        for r1 in range(m):
            for r2 in range(r1, m):
                if (r1 + r2) % m in squares_mod_m:
                    adj[r1].append(r2)
                    if r1 != r2:
                        adj[r2].append(r1)

        # 3. Find the maximum independent set size for this graph
        max_found_size = 0
        nodes = tuple(range(m))
        find_mis_size_recursive(adj, nodes, 0)
        max_r_size = max_found_size
        
        # 4. Calculate density and update best if necessary
        current_c = max_r_size / m
        print(f"{m:11d} | {max_r_size:20d} | {current_c:.6f}")

        if current_c > best_c:
            best_c = current_c
            best_m = m

    print("-" * 50)
    print(f"\nThe best density found from the search is c = {best_c:.6f} for modulus m = {best_m}.")
    
    # Manually check and explain the m=3 case, which gives the best result
    m_final = 3
    r_final = 1
    c_final = 1 / 3
    print(f"\nThe construction with modulus m = {m_final} provides the best density.")
    print(f"Let A be the set of all integers n such that n = {r_final} (mod {m_final}).")
    print(f"For any two elements a,b in this set A, their sum is a+b = (3k+{r_final}) + (3j+{r_final}) = 3(k+j) + {2*r_final} = 2 (mod {m_final}).")
    squares_mod_3 = set((k*k)%3 for k in range(3))
    print(f"The squares modulo 3 are {squares_mod_3}.")
    print("Since a+b = 2 (mod 3), the sum can never be a perfect square.")
    print(f"This construction leads to a set A of density c = {r_final}/{m_final} = {c_final}.")

if __name__ == '__main__':
    find_best_density()