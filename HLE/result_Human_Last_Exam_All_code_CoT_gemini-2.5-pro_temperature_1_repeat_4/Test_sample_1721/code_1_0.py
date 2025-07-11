import itertools

def solve():
    """
    This script searches for the largest density c = |R|/m of a set A
    such that A+A contains no square numbers, by examining modular constructions.
    """
    max_c = 0.0
    best_construction = {}

    # We limit the search to a reasonable range of moduli, as the search space grows as 2^m.
    # The optimal density is expected to be found for a small modulus.
    max_m_to_check = 16 

    print(f"Searching for the best density c = |R|/m for m from 2 to {max_m_to_check}...\n")

    for m in range(2, max_m_to_check + 1):
        # Step 1: Find all quadratic residues modulo m
        q_m = set(i*i % m for i in range(m))

        # Step 2: Find the largest subset R of {0,...,m-1} such that R+R has no quadratic residues
        best_R_for_m = set()
        max_len_R = 0

        # Iterate through all possible subsets R of {0, ..., m-1}
        # We iterate from largest possible R to smallest to find the max size faster
        for k in range(m, 0, -1):
            found_k = False
            for R_tuple in itertools.combinations(range(m), k):
                R = set(R_tuple)
                is_valid = True
                
                # Check the condition for R+R
                # itertools.product generates all pairs (r1, r2) from R, including (r,r)
                for r1, r2 in itertools.product(R, repeat=2):
                    if (r1 + r2) % m in q_m:
                        is_valid = False
                        break
                if not is_valid:
                    continue

                # Found a valid R of size k
                max_len_R = k
                best_R_for_m = R
                found_k = True
                break
            if found_k:
                break
        
        # Step 3: Calculate the density and update the maximum
        c_m = max_len_R / m if m > 0 else 0
        
        # Prepare the output string for the best set R
        r_str = str(sorted(list(best_R_for_m))) if best_R_for_m else "{}"

        print(f"m = {m:2d}: Q_{m} = {sorted(list(q_m))}")
        print(f"     Best R = {r_str}, |R| = {max_len_R}, Density |R|/m = {max_len_R}/{m} = {c_m:.4f}")
        print("-" * 30)

        if c_m > max_c:
            max_c = c_m
            best_construction = {'m': m, 'R': best_R_for_m, 'c': c_m, 'fraction': f"{max_len_R}/{m}"}

    print("\nSearch complete.")
    print(f"The largest density found is c = {best_construction.get('c', 0):.4f}")
    print(f"This corresponds to the fraction {best_construction.get('fraction', 'N/A')}.")
    m = best_construction.get('m')
    R = best_construction.get('R')
    if m is not None:
      print(f"A valid construction is the set of numbers n such that (n mod {m}) is in {R}.")


solve()
<<<1/3>>>