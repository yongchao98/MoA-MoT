import itertools

def main():
    """
    This function determines if the set Sigma is empty and computes the required value.

    It assumes the interpretation that A x A could mean element-wise product A . A,
    i.e., A+A is a subset of A*A. It then searches for a set A in Sigma.

    The search proceeds by looking for a set with the minimum possible maximum element.
    We iterate through candidate maximums 'm' in increasing order. For each 'm',
    we generate all subsets of {1, ..., m} that have 'm' as their largest element
    and check if they belong to Sigma.
    """
    
    # We will search for sets with a maximum element up to a limit.
    # If no set is found, we conclude that Sigma is empty.
    search_limit = 25 
    
    # Iterate through each possible maximum value 'm' for a set in Sigma.
    for m in range(1, search_limit + 1):
        
        # The other elements of a set with max 'm' must be smaller than 'm'.
        other_elements_pool = range(1, m)
        
        # Generate all subsets of the pool and add 'm' to form the candidate set A.
        for r in range(len(other_elements_pool) + 1):
            for combo in itertools.combinations(other_elements_pool, r):
                current_A = set(combo) | {m}
                
                # The sets {empty} and {2} are explicitly excluded from Sigma.
                # current_A is never empty. We skip {2}.
                if current_A == {2}:
                    continue
                
                # Check if the property holds for current_A.
                sum_set = {a + b for a in current_A for b in current_A}
                prod_set = {a * b for a in current_A for b in current_A}
                
                if sum_set.issubset(prod_set):
                    # If we find a set, we have found the minimum max element, since we
                    # iterate 'm' in increasing order.
                    # My analysis indicates no such set will be found.
                    # If one were found, this is how we would present the result.
                    found_set_str = ', '.join(map(str, sorted(list(current_A))))
                    print(f"min max({{{found_set_str}}}) = {m}")
                    return

    # After searching up to the limit, no set in Sigma was found.
    # This supports the conclusion that Sigma is empty.
    # As per the problem statement, we should return 0 in this case.
    # There is no set A, so no equation can be formed.
    print(0)

if __name__ == '__main__':
    main()