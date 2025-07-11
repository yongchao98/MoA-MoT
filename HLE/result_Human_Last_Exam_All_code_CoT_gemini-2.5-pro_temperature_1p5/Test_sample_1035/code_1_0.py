import math

def solve_bipartite_covering():
    """
    Calculates the minimum number of bipartite graphs needed to cover
    all edges of a complete graph K_n.
    """
    n = 35

    # The formula is k = ceil(log2(n))
    log_val = math.log2(n)
    result = math.ceil(log_val)

    print("The problem is interpreted as finding the minimum number of bipartite graphs (k) needed to cover all edges of the complete graph K_n.")
    print("The formula for this is: k = ceil(log_2(n))\n")
    print(f"For n = {n}, the calculation is as follows:")
    
    # Print the equation with the value of n
    print(f"k = ceil(log_2({n}))")
    
    # Print the equation with the intermediate value of log_2(n)
    print(f"k = ceil({log_val})")
    
    # Print the final result
    print(f"k = {result}")

solve_bipartite_covering()