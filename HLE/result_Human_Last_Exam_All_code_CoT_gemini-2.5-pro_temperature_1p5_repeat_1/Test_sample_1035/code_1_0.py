import math

def solve():
    """
    Calculates the minimum number of vertices in a family of bipartite graphs
    covering all the edges of the complete graph K_n for n=35.
    The formula is S(n) = sum_{k=1 to n} ceil(log2(k)).
    """
    n = 35
    
    terms = []
    total_sum = 0
    
    # Calculate each term from k=1 to n
    for k in range(1, n + 1):
        # log2(1) is 0, so ceil(log2(1)) is 0.
        term = math.ceil(math.log2(k)) if k > 0 else 0
        terms.append(term)
        total_sum += term
        
    # Build the equation string
    equation_str = " + ".join(map(str, terms))
    
    # Print the final equation
    print(f"The minimum number of vertices is the result of the following sum:")
    print(f"{equation_str} = {total_sum}")

solve()