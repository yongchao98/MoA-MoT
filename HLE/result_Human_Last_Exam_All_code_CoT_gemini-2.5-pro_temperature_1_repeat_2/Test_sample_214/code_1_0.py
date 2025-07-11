import sys

def solve():
    """
    This script calculates the number of open sets in the specified topological space.
    The problem reduces to counting the number of antichains in the divisibility poset on S = {1, 2, ..., 150}.
    """
    
    # Set a higher recursion limit for the recursive function.
    # The maximum depth of recursion is the size of the initial set, which is 149.
    try:
        sys.setrecursionlimit(250) 
    except (ValueError, RuntimeError): # Some environments might restrict this
        pass

    # Memoization cache to store results for subproblems, preventing re-computation.
    # Keys are sorted tuples of numbers, values are the count of antichains.
    memo = {}

    def count_antichains(s_tuple):
        """
        Recursively counts the number of antichains in a set of numbers
        ordered by divisibility. The input s_tuple is expected to be a sorted tuple.
        """
        # Base case: The empty set has one antichain, which is the empty set itself.
        if not s_tuple:
            return 1
        
        # If the result for this subproblem is already cached, return it.
        if s_tuple in memo:
            return memo[s_tuple]

        # Pick a minimal element 'm'. Since the tuple is sorted, the first element
        # is the smallest number and is always a minimal element in the divisibility poset.
        m = s_tuple[0]
        
        # Partition the antichains into two sets and count them recursively:
        # 1. Antichains that do NOT contain 'm'.
        #    These are the antichains of the set s_tuple without 'm'.
        s_without_m_tuple = s_tuple[1:]
        
        # 2. Antichains that DO contain 'm'.
        #    If an antichain contains 'm', it cannot contain any multiple of 'm'.
        #    Since 'm' is the minimum, we only need to exclude its multiples.
        s_without_neighbors_list = [k for k in s_without_m_tuple if k % m != 0]
        s_without_neighbors_tuple = tuple(s_without_neighbors_list)
        
        # Apply the recursive formula: Count(S) = Count(S - {m}) + Count(S - N(m))
        # where N(m) are elements comparable to m.
        result = count_antichains(s_without_m_tuple) + count_antichains(s_without_neighbors_tuple)
        
        # Cache the result before returning.
        memo[s_tuple] = result
        return result

    # The antichains of {1, ..., 150} are {1} union the antichains of {2, ..., 150}.
    # We calculate the number of antichains in {2, ..., 150}.
    initial_tuple = tuple(range(2, 151))
    num_antichains_without_1 = count_antichains(initial_tuple)

    # The total number of antichains is this count plus one for the antichain {1}.
    total_antichains = num_antichains_without_1 + 1

    print("The number of open sets equals the number of antichains in the poset (S, |).")
    print("This is calculated as: (Number of antichains in {2, ..., 150}) + 1")
    print(f"\nNumber of antichains in {{2, ..., 150}} = {num_antichains_without_1}")
    print(f"Adding the antichain {{1}} gives 1 more.")
    print(f"\nTotal number of open sets = {num_antichains_without_1} + 1 = {total_antichains}")

solve()