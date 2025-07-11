def solve():
    """
    Determines the number of nonzero terms, up to and including the term in x^-100,
    in the asymptotic expansion of the function f(x).
    """
    # c stores the coefficients c_n of the expansion f(x) = sum(c_n * x^-n)
    c = {}
    
    # The expansion starts from n=2, so we can define c_1 = 0.
    c[1] = 0
    
    # Calculate coefficients c_n for n from 2 to 100 using the recurrence relation.
    for n in range(2, 101):
        if n % 2 == 1:
            # For odd n, c_n = 1
            c[n] = 1
        else:
            # For even n, c_n + c_{n/2} = 1
            c[n] = 1 - c[n // 2]
            
    # Count the number of non-zero terms for odd and even n up to 100.
    odd_nonzero_count = 0
    even_nonzero_count = 0
    
    for n in range(2, 101):
        if c[n] != 0:
            if n % 2 == 1:
                odd_nonzero_count += 1
            else:
                even_nonzero_count += 1
                
    total_nonzero_count = odd_nonzero_count + even_nonzero_count
    
    # The problem asks to output the numbers in the final equation.
    # We present the sum of counts from odd and even terms.
    print(f"Number of non-zero terms for odd n (up to 100): {odd_nonzero_count}")
    print(f"Number of non-zero terms for even n (up to 100): {even_nonzero_count}")
    print(f"Total number of non-zero terms = {odd_nonzero_count} + {even_nonzero_count} = {total_nonzero_count}")

solve()