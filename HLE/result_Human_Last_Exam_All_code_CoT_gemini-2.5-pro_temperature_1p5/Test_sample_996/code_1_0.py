import math

def solve():
    """
    This function provides the determined rational numbers (a,b,c,d,e,f).

    The analysis leads to the following conclusions:
    Part A: alpha(k) is the maximum length of E(j) for j < k. The length L(k) follows a recursion
    that results in alpha(k) growing as Theta(log k). This means k is raised to the power a=0,
    log(k) is raised to the power b=1, and log(log(k)) is raised to the power c=0.
    So, (a,b,c) = (0,1,0).

    Part B: beta(k) is the count of numbers j < k with zero oddity. The count of these numbers,
    denoted beta_n for numbers up to T_n = 2^(2^n), follows the recurrence beta_{n+1} = beta_n + beta_n^2.
    This leads to an asymptotic growth of beta(k) as Theta(k^(2/3)). This implies k is raised to
    the power d=2/3, log(k) is raised to the power e=0, and log(log(k)) is raised to the power f=0.
    So, (d,e,f) = (2/3,0,0).
    """

    a = 0
    b = 1
    c = 0
    # For d, we use a string to represent the fraction 2/3 as requested.
    d_str = "2/3"
    e = 0
    f = 0
    
    # We are asked to write a,b,c,d,e,f separated by commas with no spaces.
    # The output format will have '2/3' for d.
    result = f"{a},{b},{c},{d_str},{e},{f}"
    print(result)

solve()