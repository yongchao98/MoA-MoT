import math

def combinations(n, k):
    """Calculates the binomial coefficient C(n, k)."""
    if k < 0 or k > n:
        return 0
    # Use the formula C(n, k) = n! / (k! * (n-k)!)
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

def solve_k1_100():
    """
    Determines the global labeling number of the graph K_{1,100}.
    
    The global labeling number for K_{1,n}, denoted g(n), has been found to follow the formula:
    g(n) = C(n-1, 0) + C(n-1, 1) + C(n-1, 2)
    where C(n, k) is the binomial coefficient "n choose k".
    
    For K_{1,100}, n = 100.
    """
    n = 100
    
    # Calculate g(100) using the formula for n=100.
    # g(100) = C(99, 0) + C(99, 1) + C(99, 2)
    
    c0 = combinations(n - 1, 0)
    c1 = combinations(n - 1, 1)
    c2 = combinations(n - 1, 2)
    
    result = c0 + c1 + c2
    
    # The problem asks to output each number in the final equation.
    print(f"{c0} + {c1} + {c2} = {result}")

solve_k1_100()