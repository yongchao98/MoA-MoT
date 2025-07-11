def solve():
    """
    This function determines the number of non-zero terms in the asymptotic expansion
    of f(x) up to the x^-100 term.

    The coefficients c_k of the expansion f(x) = sum(c_k * x^-k) satisfy the
    following recurrence relations:
    - c_1 = 0
    - c_2 = 1
    - c_k = 1, for odd k >= 3
    - c_k = 1 - c_{k/2}, for even k >= 4
    """
    
    # A dictionary to store the coefficients c_k
    c = {}
    
    # Initialize the first two coefficients based on the recurrence
    c[1] = 0
    c[2] = 1
    
    # Calculate coefficients from k=3 to 100
    for k in range(3, 101):
        if k % 2 != 0:  # k is odd
            c[k] = 1
        else:  # k is even
            c[k] = 1 - c[k // 2]
            
    # Count the number of non-zero terms for k from 2 to 100
    odd_nonzero_count = 0
    even_nonzero_count = 0
    
    for k in range(2, 101):
        if c[k] != 0:
            if k % 2 != 0:
                odd_nonzero_count += 1
            else:
                even_nonzero_count += 1
                
    total_nonzero_count = odd_nonzero_count + even_nonzero_count
    
    print(f"Number of non-zero terms for odd k (3 <= k <= 99): {odd_nonzero_count}")
    print(f"Number of non-zero terms for even k (2 <= k <= 100): {even_nonzero_count}")
    print(f"Total number of non-zero terms = {odd_nonzero_count} + {even_nonzero_count} = {total_nonzero_count}")

solve()