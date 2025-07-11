def solve():
    """
    Calculates the minimum number of vertices in a family of bipartite graphs
    covering all the edges of the complete graph K_n.
    """
    n = 35

    # The problem is to find the biclique cover size of K_n, which is known to be
    # the sum of integers from 2 to n.
    # The formula for this sum is n * (n + 1) / 2 - 1.

    # Calculate the sum of the first n integers
    sum_first_n = n * (n + 1) // 2

    # The result is this sum minus 1
    result = sum_first_n - 1

    # Print the explanation and the final equation with all numbers
    print("The minimum number of vertices is found using the formula: (n * (n + 1) / 2) - 1")
    print(f"For n = {n}:")
    # Show the components of the calculation
    print(f"The sum of integers from 1 to {n} is: ({n} * ({n} + 1)) / 2 = {sum_first_n}")
    print(f"The final calculation is: {sum_first_n} - 1 = {result}")
    print("\nFinal Answer:")
    print(result)

solve()