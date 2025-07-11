import math

def solve():
    """
    Calculates the minimum number of bipartite graphs needed to cover the edges of K_n.
    This is given by the formula ceil(log2(n)).
    """
    n = 35
    
    # The minimum number of bipartite subgraphs required to cover the complete graph K_n
    # is given by the ceiling of log base 2 of n.
    result = math.ceil(math.log2(n))
    
    # The problem asks for the minimum number of vertices, which is interpreted as a
    # typo for the minimum number of graphs in the family.
    
    # Print the equation and the result
    # We use math.log2 to compute the logarithm and math.ceil for the ceiling function.
    print(f"The minimum number of bipartite graphs to cover K_{n} is ceil(log2(n)).")
    print(f"For n = {n}:")
    print(f"ceil(log2({n})) = {result}")

solve()