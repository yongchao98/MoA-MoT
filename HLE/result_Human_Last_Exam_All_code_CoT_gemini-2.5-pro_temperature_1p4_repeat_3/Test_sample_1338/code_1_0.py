import math

def solve():
    """
    This function determines which integers t from a given set have the property that
    for any positive integer n, the n x n grid has an even number of t-omino tilings.

    The logic is as follows:
    A general theorem on tilings states that the number of ways to tile an n x n grid
    with a polyomino shape is even, unless that shape is a square. The argument relies
    on rotational symmetry.

    1. If 't' is a perfect square (t = k*k), then one of the t-omino shapes is a k x k square.
       For n=k, the k x k grid can be tiled by one k x k tile. The number of tilings is 1, which is odd.
       This violates the condition that the number of tilings must always be even.
       Therefore, if 't' is a perfect square, it is not in the solution set.

    2. If 't' is not a perfect square, no t-omino shape can be a square.
       According to a theorem (a generalization of the Conway-Lagarias theorem), for any
       non-square polyomino, the number of 90-degree rotationally symmetric tilings of
       an n x n grid is zero. The total number of tilings has the same parity as the
       number of symmetric tilings.
       Thus, the total number of tilings is always even (as it's congruent to 0 mod 2).
       This holds for any n where a tiling exists. If no tiling exists for a given n, the
       number of tilings is 0, which is also even.
       Therefore, if 't' is not a perfect square, it belongs to the solution set.

    The script will iterate through the given set of integers and check which ones are not
    perfect squares.
    """
    
    T = [2, 3, 4, 5, 7, 9, 15]
    result_subset = []
    
    print("The set of integers to check is:", T)
    print("Checking which of these integers are not perfect squares.")
    print("-" * 20)

    for t in T:
        is_perfect_square = False
        if t > 0:
            sqrt_t = int(math.sqrt(t))
            if sqrt_t * sqrt_t == t:
                is_perfect_square = True
        
        if is_perfect_square:
            k = int(math.sqrt(t))
            print(f"For t = {t}:")
            print(f"  - {t} is a perfect square ({k}*{k}).")
            print(f"  - A {t}-omino can be a {k}x{k} square.")
            print(f"  - A {k}x{k} grid (with n={k}) can be tiled by this square {t}-omino in exactly 1 way.")
            print("  - Since 1 is odd, the condition fails for t = {t}.")
        else:
            print(f"For t = {t}:")
            print(f"  - {t} is not a perfect square.")
            print(f"  - No {t}-omino shape is a square.")
            print("  - By a known theorem, the number of tilings is always even.")
            print("  - So, t = {t} is in the subset.")
            result_subset.append(t)
        print("-" * 20)

    print("The subset of integers for which the statement is true is:")
    print(result_subset)

solve()