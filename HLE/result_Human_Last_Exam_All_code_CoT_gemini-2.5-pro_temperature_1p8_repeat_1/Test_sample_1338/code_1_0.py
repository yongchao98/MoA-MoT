import math

def solve():
    """
    Determines for which integers t in a given set the number of t-omino tilings
    of an n x n grid is even for any positive integer n.
    """
    T = [2, 3, 4, 5, 7, 9, 15]
    result_set = []
    
    print("Analyzing the set T = {2, 3, 4, 5, 7, 9, 15}:\n")

    for t in T:
        # Test Case 1: t is a perfect square.
        sqrt_t = math.isqrt(t)
        if sqrt_t * sqrt_t == t:
            n = sqrt_t
            print(f"For t = {t}:")
            print(f"  - t is a perfect square ({t} = {n}^2).")
            print(f"  - Consider the {n}x{n} grid. Its area is {n*n}, which equals t.")
            print(f"  - A tiling requires a single {t}-omino. The only {t}-omino that can tile a {n}x{n} grid is the {n}x{n} square itself.")
            print(f"  - Thus, there is exactly 1 tiling for the {n}x{n} grid.")
            print(f"  - Since 1 is odd, t={t} is NOT in the subset.\n")
            continue

        # Test Case 2: t=2.
        if t == 2:
            print(f"For t = {t}:")
            print("  - If n is odd, n^2 is odd and not divisible by 2. Number of tilings is 0 (even).")
            print("  - If n is even, n=2k, it's a known result that the number of domino tilings of a 2kx2k grid is even.")
            print(f"  - So for any n, the number of tilings is even.")
            print(f"  - Therefore, t={t} IS in the subset.\n")
            result_set.append(t)
            continue
            
        # Test Case 3: t is odd and not a perfect square.
        if t % 2 != 0:
            print(f"For t = {t}:")
            print(f"  - t is not a perfect square and is odd.")
            print("  - The n x n grid is a centrally symmetric region.")
            print("  - The set of t-ominoes contains tiles that are not centrally symmetric (e.g., L-shapes).")
            print("  - A theorem states that if a centrally symmetric region is tiled by a set of tiles containing at least one non-centrally symmetric tile, the number of tilings is even.")
            print(f"  - Therefore, t={t} IS in the subset.\n")
            result_set.append(t)
            continue
            
    print("Summary:")
    print(f"The subset of integers for which the statement is true is: {sorted(result_set)}")

solve()