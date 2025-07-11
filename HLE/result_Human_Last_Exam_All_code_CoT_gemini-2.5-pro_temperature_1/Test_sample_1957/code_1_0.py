import math

def solve():
    """
    This function finds the minimum value of 1000m + n.
    
    The problem asks for the minimum of 1000m + n, where n and m are positive integers
    related to a polynomial matrix F. The degree of det(F) is d = n*m.
    The condition det(F)=0 is equivalent to a geometric condition on 5 points in R^3.

    Step 1: Analyze the geometric condition.
    The condition is that the 5 points (A, B, C, D, X) are either all coplanar (S1) or
    all lie on a common double cone with apex X (S2).
    If the points are coplanar with X, they lie on a degenerate cone (a plane) with apex X.
    Thus, S1 is a subset of S2, and the overall condition S = S1 U S2 simplifies to just S2.

    Step 2: Formulate the condition S2 as a polynomial equation.
    Let the apex X be the origin. The points A, B, C, D become vectors a, b, c, d.
    These four vectors lie on a cone if their "conic liftings" are linearly dependent.
    The conic lifting of a vector v=(x,y,z) is a 6D vector L(v) = (x^2, y^2, z^2, 2xy, 2xz, 2yz).
    The four vectors L(a), L(b), L(c), L(d) are linearly dependent if the 4x6 matrix M_cone
    formed by these row vectors has rank < 4 (i.e., rank <= 3).
    This is equivalent to all 4x4 minors of M_cone being zero.

    Step 3: Calculate the degree of the defining polynomial.
    The entries of M_cone are quadratic in the coordinates of the points (e.g., (A_x - X_x)^2).
    A 4x4 minor of M_cone is a polynomial of degree 4 * 2 = 8.
    For real-valued points, the condition "all minors are zero" is equivalent to "the sum of
    squares of the minors is zero".
    Let P be this sum of squares. P = m1^2 + m2^2 + ...
    The degree of P is deg(m_i^2) = 2 * deg(m_i) = 2 * 8 = 16.
    So, the degree of det(F) is d = 16.

    Step 4: Minimize 1000m + n subject to nm = 16.
    We need to find positive integers m and n that minimize the expression.
    The pairs (m, n) whose product is 16 are:
    (1, 16), (2, 8), (4, 4), (8, 2), (16, 1).
    We calculate the value for each pair:
    - m=1, n=16: 1000*1 + 16 = 1016
    - m=2, n=8:  1000*2 + 8  = 2008
    - m=4, n=4:  1000*4 + 4  = 4004
    - m=8, n=2:  1000*8 + 2  = 8002
    - m=16, n=1: 1000*16 + 1 = 16001
    The minimum value is 1016.
    """
    
    d = 16  # Total degree of det(F)
    
    min_val = float('inf')
    best_m = -1
    best_n = -1
    
    # Iterate through all integer factors of d for m
    for m in range(1, d + 1):
        if d % m == 0:
            n = d // m
            val = 1000 * m + n
            if val < min_val:
                min_val = val
                best_m = m
                best_n = n

    print(f"The degree of the polynomial det(F) is d = n * m = 16.")
    print(f"We want to minimize the expression 1000*m + n.")
    print(f"The possible integer pairs (m, n) for m*n=16 are (1,16), (2,8), (4,4), (8,2), (16,1).")
    print(f"Let's check the values:")
    print(f"m=1, n=16 => 1000*1 + 16 = 1016")
    print(f"m=2, n=8  => 1000*2 + 8  = 2008")
    print(f"m=4, n=4  => 1000*4 + 4  = 4004")
    print(f"m=8, n=2  => 1000*8 + 2  = 8002")
    print(f"m=16, n=1 => 1000*16 + 1 = 16001")
    print(f"The minimum value is {min_val}, which occurs for m={best_m} and n={best_n}.")
    
    # The final answer is the integer value
    print(min_val)

solve()