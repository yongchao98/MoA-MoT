def count_divisors(n):
    """
    Calculates tau(n), the number of divisors of a positive integer n.
    This corresponds to the number of subgroups of the cyclic group Z_n.
    """
    if not isinstance(n, int) or n < 1:
        raise ValueError("Input must be a positive integer.")
    
    count = 0
    # A number n is divisible by i if n % i == 0.
    # We check all integers from 1 up to n.
    for i in range(1, n + 1):
        if n % i == 0:
            count += 1
    return count

# Step 1: State the problem in terms of group theory.
# The number of smooth coverings corresponds to the number of non-isomorphic
# quasi-simple covering groups of S = PSL(2, p). This, in turn, is equal to
# the number of subgroups of the Schur multiplier M(S).

# Step 2: Identify the order of the Schur multiplier.
# For the simple group S = PSL(2, p) with p > 5 a prime, the Schur multiplier
# M(S) is the cyclic group of order 2.
order_of_schur_multiplier = 2

# Step 3: Calculate the number of subgroups.
# The number of subgroups of a cyclic group of order n is tau(n),
# the number of divisors of n. We calculate tau(2).
number_of_coverings = count_divisors(order_of_schur_multiplier)

# Step 4: Output the explanation, the final equation, and the result.
# The problem requires printing the numbers in the final equation.
# The final equation is: tau(2) = 2.
print("The total number of smooth coverings is determined by the number of subgroups of the Schur Multiplier of PSL(2, p).")
print(f"For a prime p > 5, the order of the Schur Multiplier is {order_of_schur_multiplier}.")
print("The number of subgroups of a cyclic group of order n is tau(n), the number of divisors of n.")
print(f"Final calculation: tau({order_of_schur_multiplier}) = {number_of_coverings}")
print(f"The total number of such smooth coverings is {number_of_coverings}.")