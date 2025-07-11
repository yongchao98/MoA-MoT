import sys

def solve():
    """
    Searches for the best modulus m and residue set B to construct a set A
    such that A+A contains no squares, maximizing the density |A|/N.
    """
    # The search for m is limited to a reasonable number, as the search for B
    # is computationally expensive (exponential in m).
    max_m = 30
    if len(sys.argv) > 1:
        try:
            max_m = int(sys.argv[1])
            print(f"Running search for m up to {max_m}.")
        except ValueError:
            print(f"Invalid argument. Using default max_m = {max_m}.")
            
    best_c = 0.0
    best_m = 0
    best_B = []

    for m in range(2, max_m + 1):
        squares_mod_m = {pow(i, 2, m) for i in range(m)}
        
        # State for the backtracking search for modulus m
        _max_b_size_for_m = 0
        _best_b_for_m = []

        # This recursive function finds the largest B for the current m
        def find_max_b_recursive(k, current_B):
            """
            k: next integer to consider for inclusion in B
            current_B: the set B being built so far
            """
            nonlocal _max_b_size_for_m, _best_b_for_m

            # Pruning: if the current path cannot lead to a better solution, stop.
            if len(current_B) + (m - k) <= _max_b_size_for_m:
                return

            if k == m:
                if len(current_B) > _max_b_size_for_m:
                    _max_b_size_for_m = len(current_B)
                    _best_b_for_m = list(current_B)
                return

            # Case 1: Don't include k in B
            find_max_b_recursive(k + 1, current_B)

            # Case 2: Try to include k in B
            # Check if adding k is valid
            is_valid_to_add = True
            # Check the sum k+k
            if (k + k) % m in squares_mod_m:
                is_valid_to_add = False
            
            # Check sums k+b for all b already in current_B
            if is_valid_to_add:
                for b_val in current_B:
                    if (k + b_val) % m in squares_mod_m:
                        is_valid_to_add = False
                        break
            
            if is_valid_to_add:
                current_B.append(k)
                find_max_b_recursive(k + 1, current_B)
                current_B.pop() # backtrack

        find_max_b_recursive(0, [])
        
        c_m = _max_b_size_for_m / m
        
        print(f"m = {m:2d}: max |B| = {_max_b_size_for_m:2d}, c = {_max_b_size_for_m}/{m} ≈ {c_m:.4f}")

        if c_m > best_c:
            best_c = c_m
            best_m = m
            best_B = _best_b_for_m

    print("\n" + "="*40)
    print("           Search Finished")
    print("="*40)
    print(f"The largest density found is c = |B|/m.")
    print(f"The best result from this search is:")
    # Outputting the numbers in the final equation: c = |B|/m
    print(f"c = {len(best_B)}/{best_m} ≈ {best_c}")
    print(f"This was found with modulus m = {best_m}.")
    print(f"The corresponding set of residues is B = {set(best_B)}.")
    print(f"\nThis suggests that the largest possible value for c is {best_c}.")
    print(f"In this case, the optimal construction is A = {{n | n mod {best_m} is in {set(best_B)}}}.")
    print(f"For any a, a' in A, a+a' mod {best_m} is in { {(b1+b2)%best_m for b1 in best_B for b2 in best_B} }.")
    print(f"The squares mod {best_m} are { {pow(i, 2, best_m) for i in range(best_m)} }.")
    print("Since the two sets of residues are disjoint, a+a' can never be a square number.")


solve()