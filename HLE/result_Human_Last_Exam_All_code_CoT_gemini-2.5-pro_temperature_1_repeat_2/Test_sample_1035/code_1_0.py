import math

def solve_bipartite_covering():
    """
    Calculates the minimum number of bipartite graphs needed to cover the
    edges of a complete graph K_n.
    """
    n = 35

    print(f"The problem is to determine the minimum number of bipartite graphs required to cover all edges of the complete graph K_n, for n = {n}.")
    print("This quantity is known as the bipartite covering number.")
    print("\nAccording to a theorem in graph theory, this number 'k' is calculated using the formula:")
    print("k = ceil(log2(chi(G)))")
    print("where chi(G) is the chromatic number of the graph G.")
    
    print(f"\nFor the complete graph K_n, the chromatic number is n.")
    print(f"In our case, G = K_{n} with n = {n}, so chi(K_{n}) = {n}.")

    print(f"\nTherefore, we need to calculate: k = ceil(log2({n}))")
    print(f"This is the same as finding the smallest integer k such that 2^k >= {n}.")

    # Find the smallest integer k such that 2^k >= n
    k = 0
    while (2**k) < n:
        k += 1

    lower_power = k - 1
    lower_bound = 2**lower_power
    upper_bound = 2**k

    print("\nLet's examine the powers of 2 around n:")
    print(f"The number in our equation is n = {n}.")
    print(f"The power of 2 just below n is 2^{lower_power} = {lower_bound}.")
    print(f"The power of 2 just at or above n is 2^{k} = {upper_bound}.")
    
    print(f"\nSince {lower_bound} < {n} <= {upper_bound}, the smallest integer k that satisfies the condition is {k}.")

    # Also show the direct log calculation
    log_value = math.log2(n)
    final_answer = math.ceil(log_value)
    
    print("\nTo confirm with logarithms:")
    print(f"The value of log2({n}) is approximately {log_value:.4f}.")
    print(f"The ceiling of {log_value:.4f} is {final_answer}.")
    
    print("\nFinal equation and answer:")
    print(f"The minimum number of bipartite graphs is k = ceil(log2({n})) = {final_answer}.")

solve_bipartite_covering()