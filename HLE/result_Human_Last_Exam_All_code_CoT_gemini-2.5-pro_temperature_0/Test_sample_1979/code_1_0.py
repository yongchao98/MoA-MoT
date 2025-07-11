import itertools

def solve():
    """
    This script aims to solve the problem by finding the set A that satisfies the given conditions.

    Problem Interpretation:
    The condition for a set A to be in Sigma is given as A+A is a subset of A x A.
    This is interpreted as a typo for A+A is a subset of A*A, where
    A*A = {a * b | a in A, b in A} is the product set of A. This interpretation is
    the most plausible one that is mathematically well-defined and non-trivial.

    The set Sigma consists of all finite, non-empty subsets of positive integers A,
    excluding the empty set and {2}, such that for all a, b in A, the sum a+b
    can be written as a product of two elements from A.

    The task is to find the minimum value of the maximum element of A, for all A in Sigma.
    If Sigma is empty, the answer is 0.

    Method:
    The script will perform a brute-force search for such a set A. It will iterate
    through possible maximum values `max_a` for the set, from 1 upwards. For each
    `max_a`, it checks all possible non-empty subsets of {1, 2, ..., max_a} that
    contain `max_a`. If a set is found that satisfies the condition, `max_a` is the
    answer, since we are searching in increasing order of `max_a`.

    Based on a mathematical proof that no such sets exist (as outlined in the thinking steps),
    this search is expected to find no solutions. If the search completes up to a
    defined limit without finding any valid set, the script will conclude that Sigma is empty
    and output 0, as per the problem statement.
    """

    def is_solution(A):
        """Checks if a set A is in Sigma."""
        # Condition: A must not be empty or {2}
        if not A or A == {2}:
            return False

        # Condition: A+A must be a subset of A*A
        sums = {a + b for a in A for b in A}
        products = {a * b for a in A for b in A}

        return sums.issubset(products)

    # We search for a solution by iterating through the possible maximum element `max_a`.
    # If a solution is found for a given `max_a`, that `max_a` is the minimum possible
    # value for max(A) because we are searching in increasing order.
    search_limit = 16  # A reasonable upper bound for the search.
    
    found_solution = False

    for max_a in range(1, search_limit + 1):
        # To construct candidate sets A with max(A) = max_a, we take all subsets
        # of {1, ..., max_a-1} and add max_a to them.
        other_elements = list(range(1, max_a))
        for i in range(len(other_elements) + 1):
            for subset_tuple in itertools.combinations(other_elements, i):
                candidate_A = set(subset_tuple)
                candidate_A.add(max_a)

                if is_solution(candidate_A):
                    # A solution is found. Since we iterate max_a upwards,
                    # this is the minimum possible max(A).
                    print(max_a)
                    found_solution = True
                    return
    
    # If the loop completes without finding any solution, Sigma is empty.
    if not found_solution:
        print(0)

solve()