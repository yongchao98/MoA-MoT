def solve():
    """
    Let C_n be the sum of squares of coefficients for a product of n terms.
    The problem asks for C_20.
    The values follow the recurrence relation C_n = 6*C_{n-1} - 2*C_{n-2}.
    Initial values:
    C_0 = 1 (for an empty product, P(x)=1, sum of squares is 1^2=1)
    C_1 = 4 (constant term of Q_0(x))
    """
    
    # We want to find C_20
    n = 20
    
    if n == 0:
        print(1)
        return
    
    # Initialize with C_0 and C_1
    c_prev = 1  # C_0
    c_curr = 4  # C_1
    
    # Iteratively compute C_n up to n=20
    for i in range(2, n + 1):
        c_next = 6 * c_curr - 2 * c_prev
        c_prev = c_curr
        c_curr = c_next
        
    print(c_curr)

solve()