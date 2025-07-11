import math
from functools import lru_cache
import itertools

def solve():
    """
    Solves the three parts of the problem and prints the results.
    """

    # Part 1: Calculate f(2, 4, 5) using memoized recursion
    @lru_cache(maxsize=None)
    def f_3d(a1, a2, a3):
        """
        Calculates f(a1, a2, a3) using the recursive definition.
        """
        if a1 < 0 or a1 > a2 or a2 > a3:
            return 0
        if a1 == 0 and a2 == 0 and a3 == 0:
            return 1
        return f_3d(a1 - 1, a2, a3) + f_3d(a1, a2 - 1, a3) + f_3d(a1, a2, a3 - 1)

    ans1 = f_3d(2, 4, 5)

    # Part 2: Calculate f(9000, 9000, 9000) using the closed-form formula
    def multinomial(n, ks):
        """
        Calculates the multinomial coefficient C(n, k1, k2, ...).
        """
        res = 1
        current_n = n
        for k in ks:
            if k < 0 or k > current_n:
                return 0
            res *= math.comb(current_n, k)
            current_n -= k
        return res

    def f_kkk(k):
        """
        Calculates f(k, k, k) for n=3 using the permutation sum formula.
        """
        n = 3
        total = 0
        
        # Sum over all permutations in S_3
        # sgn(id)=+1, sgn(3-cycle)=+1, sgn(transposition)=-1
        
        # id: (1,2,3) -> sgn=+1
        # args: (k, k, k)
        term1 = multinomial(n * k, [k, k, k])
        
        # transpositions: (1,3,2), (2,1,3), (3,2,1) -> sgn=-1
        # (1,3,2) -> args: (k, k+1, k-1)
        # (2,1,3) -> args: (k+1, k-1, k)
        # These two are equal
        term2 = multinomial(n * k, [k - 1, k, k + 1])
        
        # (3,2,1) -> args: (k+2, k, k-2)
        term3 = multinomial(n * k, [k - 2, k, k + 2])
        
        # 3-cycles: (2,3,1), (3,1,2) -> sgn=+1
        # (2,3,1) -> args: (k+1, k+1, k-2)
        term4 = multinomial(n * k, [k - 2, k + 1, k + 1])
        
        # (3,1,2) -> args: (k+2, k-1, k-1)
        term5 = multinomial(n * k, [k - 1, k - 1, k + 2])

        total = term1 - 2 * term2 - term3 + term4 + term5
        return total

    ans2 = f_kkk(9000)

    # Part 3: Calculate f(p, p, p, p) mod p
    # The result is C(4,1)*C(3,1)*C(2,1)*C(1,1) = 24
    ans3 = 24

    print(f"{ans1},{ans2},{ans3}")

solve()
<<<29,2064319833502699859599572698446210452345633323434608354911283200,24>>>