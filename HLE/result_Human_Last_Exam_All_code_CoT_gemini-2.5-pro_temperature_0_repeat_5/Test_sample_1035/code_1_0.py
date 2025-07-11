import math

def solve():
    """
    Calculates the minimum number of bipartite graphs needed to cover the edges of K_n.
    This is equivalent to finding the smallest integer k such that 2^k >= n.
    """
    n = 35

    # We are looking for the smallest integer k that satisfies the inequality 2^k >= n.
    # We can find this by checking powers of 2.
    k = 0
    power_of_2 = 1
    while power_of_2 < n:
        power_of_2 *= 2
        k += 1

    # The final "equation" is the inequality 2^k >= n.
    # We present the values around the solution to demonstrate it.
    print(f"The problem is to find the minimum integer k such that 2^k >= {n}.")
    print(f"Checking the powers of 2:")
    print(f"For k = {k-1}, the value is 2^{k-1} = {2**(k-1)}.")
    print(f"For k = {k}, the value is 2^{k} = {2**k}.")
    print(f"Since {2**(k-1)} < {n} and {2**k} >= {n}, the smallest integer k that satisfies the condition is {k}.")
    print(f"Thus, the minimum number of bipartite graphs required is {k}.")

solve()