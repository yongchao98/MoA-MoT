def solve_asymptotic_expansion():
    """
    Calculates the number of non-zero terms in the asymptotic expansion of f(x)
    up to the x^-100 term.
    """
    # The maximum power for the x^-k term
    N = 100

    # a[k] will store the coefficient of the x^-k term
    a = {}

    # From the recurrence relation a_1 = 0
    a[1] = 0

    # Calculate coefficients a_k for k from 2 to N
    for k in range(2, N + 1):
        if k % 2 == 1:
            # For odd k >= 3, the recurrence is a_k = 1
            a[k] = 1
        else:
            # For even k >= 2, the recurrence is a_k = 1 - a_{k/2}
            a[k] = 1 - a[k // 2]

    # Count the number of non-zero coefficients for odd and even k
    num_nonzero_odd = 0
    num_nonzero_even = 0
    for k in range(1, N + 1):
        if a.get(k, 0) != 0:
            if k % 2 == 1:
                num_nonzero_odd += 1
            else:
                num_nonzero_even += 1
    
    total_nonzero = num_nonzero_odd + num_nonzero_even

    # Print the final results as an equation, as requested
    print("The number of non-zero terms is calculated by summing the counts from odd and even powers.")
    print("Number of non-zero terms from odd powers (k=1 to 99): {}".format(num_nonzero_odd))
    print("Number of non-zero terms from even powers (k=2 to 100): {}".format(num_nonzero_even))
    print("\nThe final equation for the total count is:")
    print("{} + {} = {}".format(num_nonzero_odd, num_nonzero_even, total_nonzero))


solve_asymptotic_expansion()