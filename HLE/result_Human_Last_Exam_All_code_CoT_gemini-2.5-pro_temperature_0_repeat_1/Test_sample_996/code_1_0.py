def solve():
    """
    This function provides the determined rational numbers (a,b,c,d,e,f).
    
    For Part A, the analysis of the recurrence relation for the length function is complex.
    The asymptotic behavior of alpha(k) is known to be Theta(log(k)/log(log(k))).
    This corresponds to (a,b,c) = (0, 1, -1).
    
    For Part B, based on the recursive definition of oddity O(k) = O(a) + O(b),
    we can set up a recurrence for beta(k), the count of numbers with zero oddity.
    Solving this recurrence shows that beta(k) is in Theta(k^1).
    This corresponds to (d,e,f) = (1, 0, 0).
    """
    
    a = 0
    b = 1
    c = -1
    d = 1
    e = 0
    f = 0
    
    # The final answer should be printed in the specified format.
    # The problem asks for the numbers separated by commas with no spaces.
    # Example format: 0,1/2,0,1,-1,0
    # My determined numbers are all integers.
    
    print(f"{a},{b},{c},{d},{e},{f}")

solve()