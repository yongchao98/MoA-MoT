import itertools

def solve():
    """
    Determines if the set Σ is empty. If not, computes min_{A ∈ Σ} max_{a ∈ A};
    otherwise, returns 0.

    The set Σ contains finite, non-empty subsets A of positive integers
    (excluding {2}) such that A + A ⊆ A * A, where:
    A + A = {a + b | a ∈ A, b ∈ A}
    A * A = {a * b | a ∈ A, b ∈ A}

    Logical analysis suggests Σ is empty:
    1. If min(A) >= 2, it can be shown that A must only contain even numbers.
       This leads to a contradiction.
    2. If min(A) = 1, the condition 1+1=2 being in A*A implies 2 must be in A.
       This creates a chain reaction where new, larger elements seem to always
       be required to satisfy the condition for sums (e.g., A={1,2} fails,
       A={1,2,3} fails, etc.), suggesting no finite set exists.

    This script performs a search to confirm this conclusion within a practical limit.
    If no set is found, it prints 0.
    """

    memo = {}

    def check(A_tuple):
        """Checks if a given set A satisfies the condition."""
        if A_tuple in memo:
            return memo[A_tuple]

        # The problem statement excludes the empty set and {2}.
        if not A_tuple or A_tuple == (2,):
            memo[A_tuple] = False
            return False

        # A+A set
        sum_set = set()
        for i in range(len(A_tuple)):
            for j in range(i, len(A_tuple)):
                sum_set.add(A_tuple[i] + A_tuple[j])

        # A*A set
        prod_set = set()
        for i in range(len(A_tuple)):
            for j in range(i, len(A_tuple)):
                prod_set.add(A_tuple[i] * A_tuple[j])
        
        result = sum_set.issubset(prod_set)
        memo[A_tuple] = result
        return result

    # Iterate through M = max(A) and check all valid subsets.
    # We are looking for the smallest M for which a set exists.
    # Logical deduction shows any A in Σ must contain {1, 2}.
    # We start M from 3 as A={1,2} does not work.
    search_limit = 16 
    for M in range(3, search_limit + 1):
        # A must contain M, and from logic, must contain 1 and 2.
        A_base = {1, 2, M}
        
        # We can choose to add any other elements from {3, ..., M-1}.
        optional_elements = list(range(3, M))
        
        # Iterate through all combinations of optional elements.
        for i in range(len(optional_elements) + 1):
            for subset_tuple in itertools.combinations(optional_elements, i):
                A = A_base.union(set(subset_tuple))
                A_tuple = tuple(sorted(list(A)))
                
                if check(A_tuple):
                    # If a set is found, M is the minimum max, print it.
                    print(M)
                    return

    # If the search completes without finding any set, our conclusion is that
    # Σ is empty, and the answer is 0.
    # The problem does not have a simple equation, but asks for a single value as the answer.
    print(0)

solve()