import math

def solve_bipartite_covering(n):
    """
    Calculates and prints the minimum number of bipartite graphs 
    needed to cover the edges of K_n.
    """
    
    # Calculate the minimum number of bipartite graphs, k.
    # k is the smallest integer such that 2^k >= n.
    # This is equivalent to k = ceil(log2(n)).
    k = math.ceil(math.log2(n))

    print(f"For a complete graph K_n with n = {n} vertices:")
    print("We need to find the minimum number of bipartite graphs, k, to cover all its edges.")
    print(f"This number is the smallest integer k satisfying 2^k >= {n}.")
    
    # We can also find this by checking powers of 2.
    power = 0
    result = 1
    while result < n:
        power += 1
        result *= 2
    
    print("\nChecking powers of 2:")
    print(f"2^({power - 1}) = {result // 2} (which is less than {n})")
    print(f"2^{power} = {result} (which is greater than or equal to {n})")
    print(f"So, the smallest integer power required is {power}.")

    print("\nThe formal calculation is k = ceil(log2(n)).")
    log_val = math.log2(n)
    print(f"log2({n}) is approximately {log_val:.4f}.")
    
    print("\nThe final equation is:")
    print(f"ceil(log2({n})) = {k}")

# Run the solution for n=35
solve_bipartite_covering(35)