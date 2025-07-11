import math

def get_chi(a, b):
    """Computes the Euler characteristic of M(a,b) = M(a) x M(b)."""
    return 4 * (1 - a) * (1 - b)

def find_lexicographically_least_tuple():
    """
    Finds the lexicographically least tuple (a_1, b_1, ..., a_l, b_l)
    satisfying the problem's conditions by performing a systematic search.
    """
    # Generate a sorted list of candidate pairs (a,b)
    # where a, b are not 1, and a <= b.
    # The list is sorted lexicographically, which is key to finding the minimal tuple.
    max_g = 10 # A reasonable search space for genus, can be increased if needed.
    candidates = []
    for a in range(max_g):
        if a == 1:
            continue
        for b in range(a, max_g):
            if b == 1:
                continue
            candidates.append({'pair': (a, b), 'chi': get_chi(a, b)})

    # As derived in the plan, l must be odd and greater than 1. Start with l=3.
    l = 3
    while True:
        target_chi_sum = 2 * (l - 1)
        
        # Memoization for the recursive search to avoid recomputing branches.
        memo = {}

        def find_combo(target, k, start_index):
            """
            Recursively search for a multiset of k pairs from the candidates list
            (starting from start_index) whose chi values sum to the target.
            """
            if (target, k, start_index) in memo:
                return memo[(target, k, start_index)]

            if k == 0:
                return [] if target == 0 else None
            
            # Iterate through candidate pairs. The outer loop will find the first component
            # of the lexicographically smallest list of pairs.
            for i in range(start_index, len(candidates)):
                cand = candidates[i]
                new_target = target - cand['chi']
                
                # The recursive call uses `i` for start_index to allow the same pair
                # to be chosen again (forming a multiset).
                result = find_combo(new_target, k - 1, i)
                
                if result is not None:
                    # Found a valid combination.
                    solution = [cand] + result
                    memo[(target, k, start_index)] = solution
                    return solution
            
            memo[(target, k, start_index)] = None
            return None

        # Start the search for the current l.
        solution_pairs_dicts = find_combo(target_chi_sum, l, 0)
        
        if solution_pairs_dicts:
            # A solution has been found. Since we test l in increasing order, this l is minimal.
            # Since the search for pairs is on a lexicographically sorted list, the first solution
            # found will correspond to the lexicographically smallest tuple.
            
            # As requested, output the numbers in the final equation.
            chi_values = [p['chi'] for p in solution_pairs_dicts]
            equation_str = " + ".join(f"({c})" if c < 0 else str(c) for c in chi_values)
            equation_str += f" = {target_chi_sum}"
            print(equation_str)
            
            # Construct and print the final tuple as requested.
            final_tuple = []
            for p in solution_pairs_dicts:
                final_tuple.extend(p['pair'])
                
            print(str(tuple(final_tuple)).replace(" ", ""))

            return

        l += 2 # Move to the next odd number for l.

# Execute the search and print the results.
find_lexicographically_least_tuple()