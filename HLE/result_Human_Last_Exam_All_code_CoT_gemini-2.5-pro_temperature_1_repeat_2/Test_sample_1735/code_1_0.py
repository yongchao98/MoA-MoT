def solve():
    """
    This function prints the letters of the true statements in sorted order.
    A) The optimal solution can be found in O(n log n) time. (True)
    B) For any input sequence, there exists an optimal solution where all negative numbers precede all positive numbers. (False)
    C) The optimal maximum cumulative impact is invariant under negation of all elements. (True)
    D) There exists an O(n) algorithm to verify if a given permutation is optimal. (False)
    E) The problem has a unique solution if and only if all numbers are distinct and non-zero. (False)
    F) For any sequence, an optimal solution can be constructed by greedily minimizing C(k) at each step k. (False)
    G) The optimal solution remains optimal if we append any number with magnitude larger than max{C(k)}. (False)
    H) There exists a sequence where the optimal solution requires alternating positive and negative numbers. (True)
    I) The optimal solution minimizes the variance of the cumulative sums sequence. (False)
    J) For any sequence containing n pairs of numbers {x, -x}, at least n-1 such pairs must be adjacent in any optimal solution. (True)
    """
    true_statements = ['A', 'C', 'H', 'J']
    # The question asks for the letter choices of the true statements in sorted order.
    # The list is already sorted.
    print(", ".join(true_statements))

solve()