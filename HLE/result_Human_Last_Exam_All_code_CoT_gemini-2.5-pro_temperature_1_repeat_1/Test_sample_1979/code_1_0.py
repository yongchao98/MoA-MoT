import itertools

def solve_sigma_problem():
    """
    Searches for a set A that satisfies the conditions for being in Sigma,
    assuming the condition is A+A is a subset of A*A.

    The problem is to find min(max(a) for a in A) for all A in Sigma.
    Sigma is the set of all finite, non-empty subsets of positive integers A,
    such that A+A is a subset of A*A, excluding the sets empty_set and {2}.

    Our analysis shows that any such set A must contain 1 and 2.
    We search for the smallest possible maximum element, m.
    """
    # We search for m = max(A) starting from 3.
    # m=1 (A={1}) fails: A+A={2}, A*A={1}
    # m=2 (A={1,2}) fails: A+A={2,3,4}, A*A={1,2,4}
    # The case min(A)=2 only yields A={2}, which is excluded.
    # So we only need to search for sets containing {1, 2}.
    search_limit_m = 16 # Search for max(A) up to this limit.

    for m in range(3, search_limit_m):
        # A must contain 1, 2, and m.
        # The other elements must come from the set {3, 4, ..., m-1}.
        middle_elements = list(range(3, m))
        
        # Iterate through all possible subsets of the middle elements.
        for i in range(len(middle_elements) + 1):
            for subset_middle in itertools.combinations(middle_elements, i):
                # Construct the candidate set A
                A = {1, 2, m}
                A.update(subset_middle)
                
                # Check if A+A is a subset of A*A
                sum_set = {a + b for a in A for b in A}
                prod_set = {a * b for a in A for b in A}
                
                if sum_set.issubset(prod_set):
                    # Found a set A that is in Sigma.
                    # Since we are iterating m in increasing order, this is the minimum max.
                    print(f"Found a solution set A = {sorted(list(A))} with max(A) = {m}")
                    print(f"The minimum of the maximums is therefore {m}.")
                    # The problem asks for the single number value.
                    print(m)
                    return m
    
    # If the loop completes without finding any set, Sigma is likely empty.
    # The problem asks to return 0 in this case.
    # print(f"No solution found up to max(A) = {search_limit_m - 1}.")
    # print("Assuming Sigma is empty, the result is 0.")
    print(0)
    return 0

if __name__ == '__main__':
    solve_sigma_problem()
