import math

def solve():
    """
    Solves the problem by analyzing the asymptotic behavior of s_n.
    """

    # 1. Find the largest root R of r^4 - 88r - 57 = 0 using Newton's method.
    # f(r) = r^4 - 88r - 57
    # f'(r) = 4r^3 - 88
    # From an initial check, the root is between 4 and 5. Let's start with r=4.5.
    r = 4.5
    for _ in range(10): # 10 iterations are more than enough for precision
        f_r = r**4 - 88*r - 57
        fp_r = 4*r**3 - 88
        if fp_r == 0:
            break
        r = r - f_r / fp_r
    
    R = r
    
    # 2. Calculate the coefficient c1 for the dominant term in a_n.
    # a_n is approx. c1 * R^n
    # The formula for c1 is -R^3 / (4R^3 - 88).
    c1 = -R**3 / (4*R**3 - 88)

    # 3. Determine the sign of the coefficient of the dominant term in s_n.
    # The dominant term of s_n is proportional to c1*(1-3*c1) * R^(2n-1).
    # We check the sign of C = c1*(1-3*c1).
    # The denominator R^2-1 is positive since R > 4.
    C = c1 * (1 - 3 * c1)
    
    # 4. Based on the sign of C, decide the result.
    if C < 0:
        # s_n is negative for large n, so ln(s_n) is not real. Limit does not exist.
        print(-1)
    else:
        # The limit is 2 * ln(R). Calculate the required value.
        limit_val = 2 * math.log(R)
        result = 10000 * limit_val
        print(math.floor(result))

solve()