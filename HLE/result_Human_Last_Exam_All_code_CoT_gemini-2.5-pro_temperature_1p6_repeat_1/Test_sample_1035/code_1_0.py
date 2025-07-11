import math

def solve():
    """
    Calculates the minimum number of bipartite graphs needed to cover the
    edges of a complete graph K_n.
    """
    n = 35

    # The problem is to find the minimum number of bipartite graphs
    # whose edges cover the complete graph K_n.
    # The solution is given by the formula: k = ceil(log2(n)).

    # Calculate the result
    min_num_graphs = math.ceil(math.log2(n))

    # For the explanation, find the bounding powers of 2.
    log_val = math.log2(n)
    lower_power = math.floor(log_val)
    upper_power = math.ceil(log_val)
    lower_bound = 2**lower_power
    upper_bound = 2**upper_power

    print(f"The problem is to find the minimum number of bipartite graphs that cover all edges of K_n for n = {n}.")
    print("This quantity is given by the formula k = ceil(log2(n)).")
    print("\nCalculation steps:")
    print(f"1. The number of vertices is n = {n}.")
    print(f"2. We need to find the smallest integer k such that 2^k is greater than or equal to n.")
    print(f"3. We check the powers of 2 around n:")
    print(f"   2^{int(lower_power)} = {int(lower_bound)}")
    print(f"   2^{int(upper_power)} = {int(upper_bound)}")
    print(f"4. We see that {n} is between {int(lower_bound)} and {int(upper_bound)}.")
    print(f"5. The smallest integer power is {int(upper_power)}, so k must be at least {int(upper_power)}.")
    
    # Print the final equation with all numbers, as requested.
    print(f"\nThe final equation is ceil(log2({n})) = {min_num_graphs}")
    print(f"\nThus, the minimum number required is {min_num_graphs}.")

solve()