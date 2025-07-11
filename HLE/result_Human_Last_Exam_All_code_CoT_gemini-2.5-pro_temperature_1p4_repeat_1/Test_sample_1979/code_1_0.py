import itertools

def find_min_max_a():
    """
    This function searches for a set A that satisfies the given mathematical conditions.

    The problem defines a set Sigma of all finite, non-empty subsets A of positive
    integers (excluding {2}) such that A + A is a subset of A * A (product set).
    It then asks for the minimum of max(A) over all A in Sigma. If Sigma is empty,
    the answer is 0.

    Based on a mathematical proof, the only non-empty set satisfying the condition
    A+A <= A*A is A={2}. Since {2} is excluded from Sigma, Sigma is empty.
    This code performs a brute-force search for small sets to computationally verify
    this conclusion. As expected, it will find no such sets and return 0.
    """

    # We set a search limit for the maximum element in A.
    MAX_ELEMENT_LIMIT = 18

    # We iterate through possible values of max(A) (denoted as n) starting from 3,
    # because A={1}, A={2}, etc., are either invalid or excluded. The first n for which
    # a valid set A is found gives the minimum max(A).
    for n in range(3, MAX_ELEMENT_LIMIT + 1):
        
        # From the proof, we know that if a solution A exists, its minimum
        # element must be 1 or 2. We can structure our search based on this.

        # Case 1: min(A) = 1. A must contain 1 and n.
        # We choose the other elements from the range {2, 3, ..., n-1}.
        base_set_1 = {1, n}
        elements_to_choose_from_1 = list(range(2, n))
        for i in range(len(elements_to_choose_from_1) + 1):
            for subset_tuple in itertools.combinations(elements_to_choose_from_1, i):
                A = base_set_1.union(set(subset_tuple))
                
                sum_set = {a + b for a in A for b in A}
                prod_set = {a * b for a in A for b in A}
                
                if sum_set.issubset(prod_set):
                    # A solution is found. Print max(A) and terminate.
                    print(n)
                    return

        # Case 2: min(A) = 2. A must contain 2 and n.
        # We choose the other elements from the range {3, 4, ..., n-1}.
        base_set_2 = {2, n}
        elements_to_choose_from_2 = list(range(3, n))
        for i in range(len(elements_to_choose_from_2) + 1):
            for subset_tuple in itertools.combinations(elements_to_choose_from_2, i):
                A = base_set_2.union(set(subset_tuple))
                
                sum_set = {a + b for a in A for b in A}
                prod_set = {a * b for a in A for b in A}

                if sum_set.issubset(prod_set):
                    # A solution is found. Print max(A) and terminate.
                    print(n)
                    return
                    
    # If the search completes without finding any set A in Sigma, this confirms
    # (within our search limit) that Sigma is empty. The required output is 0.
    print(0)

# Execute the search function.
find_min_max_a()