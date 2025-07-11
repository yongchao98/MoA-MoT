def solve():
    """
    Calculates the number of nonzero terms up to and including x^-100 in the
    asymptotic expansion of the function f(x).
    """
    # a_coeffs will store the calculated coefficients a_n.
    a_coeffs = {}
    nonzero_count = 0

    # Loop through n from 2 to 100 to calculate each a_n.
    for n in range(2, 101):
        # Apply the recurrence relation for a_n.
        if n == 2:
            # Base case for n=2
            a_n = 1
        elif n % 2 != 0:
            # Case for odd n > 2
            a_n = 1
        else: # Case for even n > 2
            a_n = 1 - a_coeffs[n // 2]
        
        # Store the calculated coefficient.
        a_coeffs[n] = a_n

        # Increment the counter if the coefficient is not zero.
        if a_n != 0:
            nonzero_count += 1
            
    print(f"The number of non-zero terms up to x^-100 is: {nonzero_count}")

solve()