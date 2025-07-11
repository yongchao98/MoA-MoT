def get_max_prime_implicants(n):
    """
    Returns the maximal number of prime implicants of a Boolean function of n variables.
    
    This function uses pre-computed values from established literature, as there is
    no simple closed-form formula for this sequence.
    The sequence a(n) starts with {1, 2, 4, 10, 26, 76, ...} for n = 0, 1, 2, 3, 4, 5.
    """
    known_values = {
        0: 1,
        1: 2,
        2: 4,
        3: 10,
        4: 26,
        5: 76
    }
    
    if n in known_values:
        return known_values[n]
    else:
        raise NotImplementedError(f"The value for a({n}) is not available in the pre-computed list.")

# The user wants to find the value for n=4.
n = 4
result = get_max_prime_implicants(n)

# Print the final equation with the numbers involved.
print(f"a({n}) = {result}")
