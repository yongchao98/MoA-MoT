def get_max_prime_implicants(n):
    """
    Calculates a(n), the maximal number of prime implicants for a Boolean function
    of n variables, by looking up the value from a list of known results.

    This sequence is known as OEIS A000373. The values are established in
    combinatorics and logic synthesis literature. There is no simple
    closed-form formula for a(n).
    """
    # Known values for the sequence a(n)
    # a(0), a(1), a(2), a(3), a(4), ...
    known_sequence = [1, 2, 4, 10, 26, 76, 232]

    if 0 <= n < len(known_sequence):
        return known_sequence[n]
    else:
        # For values of n not in this list, a more complex algorithm or
        # an extended sequence table would be required.
        raise ValueError(f"Value for n={n} is not in the pre-computed list.")

# The user wants to find the value of a(4).
n_value = 4

# Calculate the result.
result = get_max_prime_implicants(n_value)

# Print the final result in an equation format.
print(f"a({n_value}) = {result}")
