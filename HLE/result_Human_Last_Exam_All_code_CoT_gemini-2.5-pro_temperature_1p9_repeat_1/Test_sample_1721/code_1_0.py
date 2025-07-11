import math

# Store the details of the best result found
best_c = 0
best_m = 0
best_R = []

# This is a global variable to be updated by the recursive search
# It will store the largest compatible set R found for a given m
_m_best_R = []

def find_largest_R_recursive(m, S_m, current_R, candidates):
    """
    A recursive backtracking search to find the largest set R such that
    R+R has no elements in S_m.
    - m: the modulus
    - S_m: the set of quadratic residues mod m
    - current_R: the compatible set being built in the current recursion path
    - candidates: the list of remaining candidates to try adding to current_R
    """
    global _m_best_R
    
    # Pruning: if the current set plus all remaining candidates can't beat
    # the best found so far, we can stop searching this path.
    if len(current_R) + len(candidates) <= len(_m_best_R):
        return

    # If we are at a leaf of the search (no more candidates)
    # or just found a better solution, update the best one.
    if len(current_R) > len(_m_best_R):
        _m_best_R = list(current_R)

    # Iterate through remaining candidates
    for i in range(len(candidates)):
        node = candidates[i]
        
        # Check if this new node is compatible with the current_R set
        is_compatible = True
        for member in current_R:
            if (member + node) % m in S_m:
                is_compatible = False
                break
        
        if is_compatible:
            # If compatible, recurse deeper with the new node added
            new_candidates = candidates[i+1:]
            find_largest_R_recursive(m, S_m, current_R + [node], new_candidates)


def solve_for_modulus(m):
    """
    For a given modulus m, finds the largest c = |R|/m.
    """
    global _m_best_R
    
    # Step 1: Calculate the quadratic residues S_m
    S_m = {x*x % m for x in range(m)}
    
    # Step 2: Pre-filter candidates. An element 'r' can only be in R if
    # 2*r = r+r is not a quadratic residue.
    initial_candidates = [r for r in range(m) if (2*r) % m not in S_m]
    
    # Step 3: Use backtracking to find the largest compatible set R.
    # We reset the global result for this 'm' before starting the search.
    _m_best_R = []
    find_largest_R_recursive(m, S_m, [], initial_candidates)
    
    k = len(_m_best_R)
    c = k / m if m > 0 else 0
    
    print(f"m = {m:2d}: k = {k:2d}, c = {k}/{m} ≈ {c:.4f}. R = {sorted(_m_best_R)}")
    
    return c, k, sorted(_m_best_R)

def main():
    """
    Main function to drive the search for the best c.
    """
    global best_c, best_m, best_R

    # We search moduli up to 40. This is sufficient to find the optimal value
    # that is known in the literature.
    m_max = 40
    print(f"Searching for the largest c = k/m up to m = {m_max}...\n")
    
    for m in range(3, m_max + 1):
        c, k, R = solve_for_modulus(m)
        if c > best_c:
            best_c = c
            best_m = m
            best_R = R
    
    print("\n--- Search Complete ---\n")
    print("The largest density c found is c = k/m.")
    print(f"Equation: c = {len(best_R)}/{best_m}")
    print(f"Decimal value: c ≈ {best_c}")
    print(f"This was achieved with modulus m = {best_m} and residue set R = {best_R}.")
    
    print("\n### Analysis of the Best Result ###")
    if best_m % 3 == 0 and all(r % 3 == 1 for r in best_R):
        print("The optimal construction found is A = {n | n ≡ 1 (mod 3)}.")
        print("For any two elements a1, a2 from this set, a1 = 3k+1 and a2 = 3j+1.")
        print("Their sum is a1 + a2 = 3(k+j) + 2.")
        print("A number that is 2 (mod 3) can never be a perfect square, as k^2 mod 3 is always 0 or 1.")
        print("Therefore, the set A+A contains no squares.")
        print("The density of this set is 1/3.")
    else:
        # Fallback explanation if a different result emerges
        print(f"With the set A = {{n | n mod {best_m} is in {best_R}}},")
        print("the sumset A+A contains no perfect squares.")
        print(f"The density of this set is {len(best_R)}/{best_m}.")

if __name__ == '__main__':
    main()